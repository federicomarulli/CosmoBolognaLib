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
!
! file to be inserted in pix_tools.F90 with
!  #include "pixel_routines.F90"
!
! the following routines should be compiled twice:
! once without DOI8B, and once with DOI8B
! in order to respectively support 4-byte and 8-byte pixel indexing
!
!   s pix2ang_ring
!   s pix2vec_ring
!   s ang2pix_ring
!   s vec2pix_ring
!   s pix2ang_nest
!   s pix2vec_nest
!   s ang2pix_nest
!   s vec2pix_nest
!   s nest2ring
!   s ring2nest
!   f npix2nside
!   s in_ring
!   s getdisc_ring
!   s query_disc
!   s query_strip
!   s query_polygon
!   s query_triangle
!   s xy2pix_nest
!   s pix2xy_nest
!   s neighbours_nest
!   s next_in_line_nest
!   s template_pixel_ring
!   s same_shape_pixels_ring
!   s template_pixel_nest
!   s same_shape_pixels_nest
!   s uniq2nest
!   s nest2uniq
!
! 2012-03-02: make sure that both arguments of iand and modulo 
!      are of the same type (for XL Fortran compiler)
!=======================================================================
!     pix2ang_ring
!
!     renders theta and phi coordinates of the nominal pixel center
!     for the pixel number ipix (RING scheme)
!     given the map resolution parameter nside
!
!=======================================================================
#ifdef DOI8B
  subroutine pix2ang_ring_8(nside, ipix, theta, phi)
    integer(i4b), parameter :: MKD = I8B
#else
  subroutine pix2ang_ring  (nside, ipix, theta, phi)
    integer(i4b), parameter :: MKD = I4B
#endif
    INTEGER(KIND=I4B), INTENT(IN)  :: nside
    INTEGER(KIND=MKD), INTENT(IN)  :: ipix
    REAL(KIND=DP),     INTENT(OUT) :: theta, phi

    INTEGER(KIND=I4B) ::  nl2, nl4, iring, iphi
    INTEGER(KIND=MKD) ::  npix, ncap, ip
    REAL(KIND=DP) ::  fodd, dnside
    real(kind=dp), parameter :: half = 0.500000000000000_dp
    !real(kind=dp), parameter :: one  = 1.000000000000000_dp
    !real(kind=dp), parameter :: three = 3.00000000000000_dp
    real(kind=dp), parameter :: threehalf = 1.50000000000000_dp
    character(len=*), parameter :: code = "pix2ang_ring"
    !-----------------------------------------------------------------------
#ifdef DOI8B
    if (nside <= ns_max4) then ! use faster 32-bit routine whenever possible
       call pix2ang_ring(nside, int(ipix,kind=i4b), theta, phi)
       return
    endif
#else
    if (nside > ns_max4) call fatal_error(code//"> nside out of range")
#endif
    npix = nside2npix(nside)       ! total number of points
    if (ipix <0 .or. ipix>npix-1) call fatal_error (code//"> ipix out of range")

    nl2  = 2*nside
    ncap = nl2*(nside-1_MKD) ! points in each polar cap, =0 for nside =1
    dnside = real(nside, kind=dp)

    if (ipix < ncap) then ! North Polar cap -------------

!       iring = nint( sqrt( (ipix+1) * half ), kind=MKD) ! counted from North pole
       iring = (cheap_isqrt(2*ipix+2) + 1)/2
       iphi  = ipix - 2*iring*(iring - 1_MKD)
! #ifdef DOI8B
!        ! fix round-off error appearing at large Nside
!        call correct_ring_phi(1, iring, iphi)
! #endif

!       theta = ACOS( one - (iring/dnside)**2 / three )
       theta = 2.0_dp * asin(iring / (sqrt(6.0_dp)*dnside))
       phi   = (real(iphi,kind=dp) + half) * HALFPI/iring

    elseif (ipix < npix-ncap) then ! Equatorial region ------

       ip    = ipix - ncap
       nl4   = 4*nside
       iring = INT( ip / nl4 ) + nside ! counted from North pole
       iphi  = iand(ip, nl4-1_MKD)

       fodd  = half * ( iand(iring+nside+1,1) )  ! 0 if iring+nside is odd, 1/2 otherwise
       theta = ACOS( (nl2 - iring) / (threehalf*dnside) )
       phi   = (real(iphi,kind=dp) + fodd) * HALFPI / dnside

    else ! South Polar cap -----------------------------------

       ip    = npix - ipix
!       iring = nint( sqrt( ip * half ), kind=MKD)     ! counted from South pole
       iring = (cheap_isqrt(2*ip) + 1) / 2
       iphi  = 2*iring*(iring + 1_MKD) - ip
! #ifdef DOI8B
!        ! fix round-off error appearing at large Nside
!        call correct_ring_phi(-1, iring, iphi)
! #endif

!       theta = ACOS( (iring/dnside)**2 / three  - one)
       theta = PI - 2.d0 * asin(iring / (sqrt(6.0_dp)*dnside))
       phi   = (real(iphi,kind=dp) + half) * HALFPI/iring

    endif

    return
#ifdef DOI8B
  end subroutine pix2ang_ring_8
#else
  end subroutine pix2ang_ring
#endif

!=======================================================================
!     pix2vec_ring
!
!     renders vector (x,y,z) coordinates of the nominal pixel center
!     for the pixel number ipix (RING scheme)
!     given the map resolution parameter nside
!     also returns the (x,y,z) position of the 4 pixel vertices (=corners)
!     in the order N,W,S,E
!=======================================================================

#ifdef DOI8B
  subroutine pix2vec_ring_8(nside, ipix, vector, vertex)
    integer(i4b), parameter :: MKD = I8B
#else
  subroutine pix2vec_ring  (nside, ipix, vector, vertex)
    integer(i4b), parameter :: MKD = I4B
#endif
    INTEGER(KIND=I4B), INTENT(IN)                             :: nside
    INTEGER(KIND=MKD), INTENT(IN)                             :: ipix
    REAL(KIND=DP),     INTENT(OUT),dimension(1:)              :: vector
    REAL(KIND=DP),     INTENT(OUT),dimension(1:,1:), optional :: vertex

    INTEGER(KIND=I4B) :: nl2, nl4, iring, iphi
    INTEGER(KIND=MKD) ::  npix, ncap, ip
    REAL(KIND=DP) ::  fact1, fact2, fodd, z, sth, phi
    real(kind=dp), parameter :: half = 0.500000000000000_dp

    real(kind=DP) :: phi_nv, phi_wv, phi_sv, phi_ev, sin_phi, cos_phi
    real(kind=DP) :: z_nv, z_sv, sth_nv, sth_sv
    real(kind=DP) :: hdelta_phi
    integer(kind=I4B) :: iphi_mod, iphi_rat
    logical(kind=LGT) :: do_vertex
    integer(kind=i4b) :: diff_phi
    character(len=*), parameter :: code = "pix2vec_ring"
    !-----------------------------------------------------------------------
#ifdef DOI8B
    if (nside <= ns_max4) then ! use faster 32-bit routine whenever possible
       if (present(vertex)) then
          call pix2vec_ring(nside, int(ipix, kind=i4b), vector, vertex)
       else
          call pix2vec_ring(nside, int(ipix, kind=i4b), vector)
       endif
       return
    endif
#else
    if (nside > ns_max4) call fatal_error(code//"> nside out of range")
#endif
    npix = nside2npix(nside)       ! total number of points
    if (ipix <0 .or. ipix>npix-1) call fatal_error(code//"> ipix out of range")

    nl2   = 2*nside
    ncap  = nl2*(nside-1_MKD) ! points in each polar cap, =0 for nside =1

    do_vertex = .false.
    if (present(vertex)) then
       if (size(vertex,dim=1) >= 3 .and. size(vertex,dim=2) >= 4) then
          do_vertex = .true.
       else
          call fatal_error(code//">  vertex array has wrong size ")
       endif
    endif

    if (ipix < ncap) then ! North Polar cap -------------

!        iring = nint( sqrt( (ipix+1) * half ), kind=MKD) ! counted from North pole
       iring = (cheap_isqrt(2*ipix+2) + 1) / 2
       iphi  = ipix - 2*iring*(iring - 1_MKD)
! #ifdef DOI8B
!        ! fix round-off error appearing at large Nside
!        call correct_ring_phi(1, iring, iphi)
! #endif

       fact2 = (3.00000_dp*nside)*nside
       z =  1.0_dp - (iring / fact2) * iring
       sth = iring * sqrt( (1.d0+z)/fact2 )
       phi   = (real(iphi,kind=dp) + half) * HALFPI/iring

       if (do_vertex) then
          hdelta_phi = PI/(4.0_dp*iring)   ! half pixel width
          z_nv = 1.0_dp - (iring-1) / fact2 * (iring-1)
          z_sv = 1.0_dp - (iring+1) / fact2 * (iring+1)
          iphi_mod = MODULO(iphi, iring) ! in {0,1,... iring-1}
          iphi_rat = (iphi) / iring      ! in {0,1,2,3}
          phi_nv = 0.0_dp
          if (iring > 1) phi_nv = HALFPI * (iphi_rat +  iphi_mod   /real(iring-1,kind=dp))
          phi_sv                = HALFPI * (iphi_rat + (iphi_mod+1)/real(iring+1,kind=dp))
          diff_phi = 3 ! both phi_nv and phi_sv different from phi
       endif


    elseif (ipix < npix - ncap) then ! Equatorial region ------

       nl4 = 4*nside
       ip    = ipix - ncap
       iring = INT( ip / nl4 ) + nside ! counted from North pole
       iphi  = iand(ip, nl4-1_MKD)

       fact1 =  1.50000_dp*nside
       fodd  = half * ( iand(iring+nside+1,1) )  ! 0 if iring+nside is odd, 1/2 otherwise
       z = (nl2 - iring) / fact1
       sth = SQRT((1.0_dp-z)*(1.0_dp+z))
       phi   = (real(iphi,kind=dp) + fodd) * HALFPI / nside

       if (do_vertex) then
          fact2 = (3.00000_dp*nside)*nside
          hdelta_phi = PI/(4.0_dp*nside)   ! half pixel width
          phi_nv = phi
          phi_sv = phi
          z_nv = (nl2 - iring +1) / fact1
          z_sv = (nl2 - iring -1) / fact1
          diff_phi = 0 ! phi_nv = phi_sv = phi
          if (iring == nside) then ! northern transition
             z_nv = 1.0_dp - (nside-1) / fact2 * (nside-1)
             iphi_mod = iand(iphi, nside-1) ! in {0,1,... nside-1}
             iphi_rat = (iphi) / nside      ! in {0,1,2,3}
             if (nside > 1) then
                phi_nv = HALFPI * (iphi_rat +  iphi_mod /real(nside-1,kind=dp))
                diff_phi = 1
             endif
          elseif (iring == 3*nside) then ! southern transition
             z_sv = -1.0_dp + (nside-1) / fact2 * (nside-1)
             iphi_mod = iand(iphi, nside-1) ! in {0,1,... iring-1}
             iphi_rat = (iphi) / nside      ! in {0,1,2,3}
             if (nside > 1) then
                phi_sv = HALFPI * (iphi_rat +  iphi_mod   /real(nside-1,kind=dp))
                diff_phi = 2
             endif
          endif
       endif

    else ! South Polar cap -----------------------------------

       ip    = npix - ipix
!        iring = nint( sqrt( ip * half ), kind=MKD)     ! counted from South pole
       iring = (cheap_isqrt(2*ip) + 1)/2
       iphi  = 2*iring*(iring + 1_MKD) - ip
! #ifdef DOI8B
!        ! fix round-off error appearing at large Nside
!        call correct_ring_phi(-1, iring, iphi)
! #endif

       fact2 = (3.00000_dp*nside)*nside
       z = -1.0_dp + (iring / fact2) * iring
       sth = iring * sqrt( (1.0_dp - z)/fact2 )
       phi   = (real(iphi,kind=dp) + half) * HALFPI/iring

       if (do_vertex) then
          hdelta_phi = PI/(4.0_dp*iring)   ! half pixel width
          z_nv = -1.0_dp + (iring+1)/ fact2 * (iring+1)
          z_sv = -1.0_dp + (iring-1)/ fact2 * (iring-1)
          iphi_mod = MODULO(iphi, iring) ! in {0,1,... iring-1}
          iphi_rat = (iphi) / iring      ! in {0,1,2,3}
          phi_nv                = HALFPI * (iphi_rat + (iphi_mod+1)/real(iring+1,kind=dp))
          phi_sv = 0.0_dp
          if (iring > 1) phi_sv = HALFPI * (iphi_rat +  iphi_mod   /real(iring-1,kind=dp))
          diff_phi = 3
       endif

    endif

    ! pixel center
    !!!!sth = SQRT((1.0_dp-z)*(1.0_dp+z))
    cos_phi = cos(phi)
    sin_phi = sin(phi)
    vector(1) = sth * cos_phi
    vector(2) = sth * sin_phi
    vector(3) = z

    if (do_vertex) then
       ! west vertex
       phi_wv      = phi - hdelta_phi
       vertex(1,2) = sth * COS(phi_wv)
       vertex(2,2) = sth * SIN(phi_wv)
       vertex(3,2) = z

       ! east vertex
       phi_ev      = phi + hdelta_phi
       vertex(1,4) = sth * COS(phi_ev)
       vertex(2,4) = sth * SIN(phi_ev)
       vertex(3,4) = z

       ! north and south vertices
       sth_nv = SQRT((1.0_dp-z_nv)*(1.0_dp+z_nv))
       sth_sv = SQRT((1.0_dp-z_sv)*(1.0_dp+z_sv))
       if (diff_phi == 0) then
          vertex(1,1) = sth_nv * cos_phi
          vertex(2,1) = sth_nv * sin_phi
          vertex(1,3) = sth_sv * cos_phi
          vertex(2,3) = sth_sv * sin_phi
       else
          vertex(1,1) = sth_nv * COS(phi_nv)
          vertex(2,1) = sth_nv * SIN(phi_nv)
          vertex(1,3) = sth_sv * COS(phi_sv)
          vertex(2,3) = sth_sv * SIN(phi_sv)
       endif
       vertex(3,1) = z_nv
       vertex(3,3) = z_sv
    endif

    return
#ifdef DOI8B
  end subroutine pix2vec_ring_8
#else
  end subroutine pix2vec_ring
#endif

!=======================================================================
!     ang2pix_ring
!
!     renders the pixel number ipix (RING scheme) for a pixel which contains
!     a point on a sphere at coordinates theta and phi, given the map
!     resolution parameter nside
!=======================================================================
#ifdef DOI8B
  subroutine ang2pix_ring_8(nside, theta, phi, ipix)
    integer(i4b), parameter :: MKD = I8B
    integer(i4b), parameter :: mykind = I4B ! replace with i8b to reach nside=2^29
#else
  subroutine ang2pix_ring  (nside, theta, phi, ipix)
    integer(i4b), parameter :: MKD = I4B
    integer(i4b), parameter :: mykind = I4B
#endif
    INTEGER(KIND=I4B), INTENT(IN)  :: nside
    REAL(KIND=DP),     INTENT(IN)  :: theta, phi
    INTEGER(KIND=MKD), INTENT(OUT) :: ipix

    INTEGER(mykind)   ::  nl4, jp, jm, ir, ip
    REAL(KIND=DP)     ::  z, za, tt, tp, tmp, temp1, temp2
    INTEGER(KIND=I4B) ::  kshift, ipix4
    character(len=*), parameter :: code = "ang2pix_ring"

    !-----------------------------------------------------------------------
#ifdef DOI8B
    if (nside <= ns_max4) then ! use faster 32-bit routine whenever possible
       call ang2pix_ring(nside, theta, phi, ipix4)
       ipix = int(ipix4, kind=MKD)
       return
    endif
    if (nside <1 .or. nside > ns_max ) call fatal_error(code//"> nside out of range")
#else
    if (nside <1 .or. nside > ns_max4) call fatal_error(code//"> nside out of range")
#endif
    if (theta<0.0_dp .or. theta>pi)  then
       print*,code//"> theta : ",theta," is out of range [0, Pi]"
       call fatal_error
    endif

    z = COS(theta)
    za = ABS(z)
    tt = MODULO( phi, twopi) / halfpi  ! in [0,4)

    nl4 = 4_mykind*nside
    if ( za <= twothird ) then ! Equatorial region ------------------
       temp1 = nside*(0.50000_dp + tt)
       temp2 = nside* 0.75000_dp * z
       jp = int(temp1-temp2, kind=mykind) ! index of  ascending edge line
       jm = int(temp1+temp2, kind=mykind) ! index of descending edge line

       ir = nside + jp - jm ! in {0,2n} (ring number counted from z=2/3)
       kshift = iand(ir, 1) ! kshift=1 if ir is odd, 0 otherwise

       ip = INT( ( jp+jm - nside + kshift + 1 ) / 2 , kind=mykind) ! in {0,4n-1}
       if (ip >= nl4) ip = ip - nl4

       ipix = 2*nside*(nside-1_MKD) + nl4*int(ir,MKD) + ip

    else ! North & South polar caps -----------------------------

       tp = tt - INT(tt)      !MODULO(tt,1.0_dp)
       !tmp = nside * SQRT( 3.0_dp*(1.0_dp - za) )
       if (z > 0.0_dp) then
          tmp = sqrt(6.0_dp) * sin( theta * 0.5_dp) * nside
       else
          tmp = sqrt(6.0_dp) * cos( theta * 0.5_dp) * nside
       endif

       jp = INT(  tp          * tmp, kind=mykind) ! increasing edge line index
       jm = INT((1.0_dp - tp) * tmp, kind=mykind) ! decreasing edge line index

       ir = jp + jm + 1        ! ring number counted from the closest pole
       ip = INT( tt * ir, kind=mykind)     ! in {0,4*ir-1}
       if (ip >= 4*ir) ip = ip - 4*ir

       if (z>0._dp) then
          ipix = 2*ir*(ir-1_MKD) + ip
       else
          ipix = 3*nside*int(nl4,MKD) - 2*ir*(ir+1_MKD) + ip
       endif

    endif

    return
#ifdef DOI8B
  end subroutine ang2pix_ring_8
#else
  end subroutine ang2pix_ring
#endif
!=======================================================================
!     vec2pix_ring
!
!     renders the pixel number ipix (RING scheme) for a pixel which contains
!     a point on a sphere at coordinate vector (=x,y,z), given the map
!     resolution parameter nside
!=======================================================================
#ifdef DOI8B
  subroutine vec2pix_ring_8(nside, vector, ipix)
    integer(i4b), parameter :: MKD = I8B
    integer(i4b), parameter :: mykind = I4B ! replace with i8b to reach nside=2^29
#else
  subroutine vec2pix_ring  (nside, vector, ipix)
    integer(i4b), parameter :: MKD = I4B
    integer(i4b), parameter :: mykind = I4B
#endif
    INTEGER(KIND=I4B), INTENT(IN)                :: nside
    REAL(KIND=DP),     INTENT(IN), dimension(1:) :: vector
    INTEGER(KIND=MKD), INTENT(OUT)               :: ipix

    INTEGER(mykind)   :: nl4, jp, jm, ir, ip
    REAL(KIND=DP)     :: z, za, tt, tp, tmp, dnorm, phi,temp1,temp2
    INTEGER(KIND=I4B) :: kshift, ipix4
    character(len=*), parameter :: code = "vec2pix_ring"

    !-----------------------------------------------------------------------
#ifdef DOI8B
    if (nside <= ns_max4) then ! use faster 32-bit routine whenever possible
       call vec2pix_ring(nside, vector, ipix4)
       ipix = int(ipix4, kind=MKD)
       return
    endif
    if (nside<1 .or. nside>ns_max ) call fatal_error(code//"> nside out of range")
#else
    if (nside<1 .or. nside>ns_max4) call fatal_error(code//"> nside out of range")
#endif
    dnorm = SQRT(vector(1)**2+vector(2)**2+vector(3)**2)
    z  = vector(3) / dnorm  ! cos(theta)
    phi = 0.0_dp
    if (vector(1) /= 0.0_dp .or. vector(2) /= 0.0_dp) &
         &     phi = ATAN2(vector(2),vector(1)) ! phi in ]-pi,pi]

    za = ABS(z)
    if (phi < 0.0_dp)     phi = phi + twopi ! phi in [0,2pi[
    tt = phi / halfpi   ! in [0,4)

    nl4 = 4_mykind*nside
    if ( za <= twothird ) then ! Equatorial region ------------------
       temp1 = nside*(0.50000_dp + tt)
       temp2 = nside* 0.75000_dp * z
       jp = int(temp1-temp2, kind=mykind) ! index of  ascending edge line
       jm = int(temp1+temp2, kind=mykind) ! index of descending edge line
!        jp = INT(nside*(0.5_dp + tt - z*0.75_dp), kind=mykind) index of  ascending edge line
!        jm = INT(nside*(0.5_dp + tt + z*0.75_dp), kind=mykind) index of descending edge line

       ir = nside + jp - jm ! in {0,2n} (ring number counted from z=2/3)
       kshift = iand(ir, 1) ! kshift=1 if ir is odd, 0 otherwise

       ip = INT( ( jp+jm - nside + kshift + 1 ) / 2, kind=mykind) ! in {0,4n-1}
       if (ip >= nl4) ip = ip - nl4

       ipix = 2*nside*(nside-1_MKD) + nl4*int(ir,MKD) + ip

    else ! North & South polar caps -----------------------------

       tp = tt - INT(tt)      !MODULO(tt,1.0_dp)
       !tmp = nside * SQRT( 3.0_dp*(1.0_dp - za) )
       tmp = SQRT(vector(1)**2+vector(2)**2) / dnorm ! sin(theta)
       tmp = nside * tmp * SQRT( 3.0_dp / (1.0_dp + za) ) !more accurate

       jp = INT(  tp          * tmp, kind=mykind) ! increasing edge line index
       jm = INT((1.0_dp - tp) * tmp, kind=mykind) ! decreasing edge line index

       ir = jp + jm + 1        ! ring number counted from the closest pole
       ip = INT( tt * ir, kind=mykind)     ! in {0,4*ir-1}
       if (ip >= 4*ir) ip = ip - 4*ir

       if (z>0._dp) then
          ipix = 2*ir*(ir-1_MKD) + ip
       else
          ipix = 3*nside*int(nl4,MKD) - 2*ir*(ir+1_MKD) + ip
       endif

    endif

    return
#ifdef DOI8B
  end subroutine vec2pix_ring_8
#else
  end subroutine vec2pix_ring
#endif
!=======================================================================
!     pix2ang_nest
!
!     renders theta and phi coordinates of the nominal pixel center
!     for the pixel number ipix (NESTED scheme)
!     given the map resolution parameter nside
!=======================================================================
#ifdef DOI8B
  subroutine pix2ang_nest_8(nside, ipix, theta, phi)
    integer(i4b), parameter :: MKD = I8B
#else
  subroutine pix2ang_nest  (nside, ipix, theta, phi)
    integer(i4b), parameter :: MKD = I4B
#endif
    INTEGER(KIND=I4B), INTENT(IN)  :: nside
    INTEGER(KIND=MKD), INTENT(IN)  :: ipix
    REAL(KIND=DP),     INTENT(OUT) :: theta, phi

    INTEGER(KIND=MKD) :: npix, npface, ipf
    INTEGER(KIND=I4B) :: ip_low, ip_trunc, ip_med, ip_hi, &
         &     jrt, jr, nr, jpt, jp, kshift, nl4, scale, i, ismax
    INTEGER(KIND=I4B) :: ix, iy, face_num
    REAL(KIND=DP)     :: z, fn, fact1, fact2

    ! coordinate of the lowest corner of each face
    INTEGER(KIND=I4B), dimension(1:12) :: jrll = (/ 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 /) ! in unit of nside
    INTEGER(KIND=I4B), dimension(1:12) :: jpll = (/ 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 /) ! in unit of nside/2
    character(len=*), parameter :: code = "pix2ang_nest"
    !-----------------------------------------------------------------------
#ifdef DOI8B
    if (nside <= ns_max4) then ! use faster 32-bit routine whenever possible
       call pix2ang_nest(nside, int(ipix,kind=i4b), theta, phi)
       return
    endif
#else
    if (nside > ns_max4) call fatal_error(code//"> nside out of range")
#endif
    npix = nside2npix(nside)       ! total number of points
    if (ipix <0 .or. ipix>npix-1) call fatal_error ("ipix out of range")

    !     initiates the array for the pixel number -> (x,y) mapping
    if (pix2x(1023) <= 0) call mk_pix2xy()

    npface = nside * int(nside, kind=MKD)
    nl4    = 4_MKD*nside

    !     finds the face, and the number in the face
    face_num = ipix/npface  ! face number in {0,11}
    ipf = MODULO(ipix,npface)  ! pixel number in the face {0,npface-1}

    fn = real(nside, kind=dp)
    !fact1 = 1.0_dp/(3.0_dp*fn*fn)

    !     finds the x,y on the face (starting from the lowest corner)
    !     from the pixel number
    if (nside <= ns_max4) then
       ip_low = iand(ipf,1023_MKD)       ! content of the last 10 bits
       ip_trunc =    ipf/1024        ! truncation of the last 10 bits
       ip_med = iand(ip_trunc,1023)  ! content of the next 10 bits
       ip_hi  =      ip_trunc/1024   ! content of the high weight 10 bits

       ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
       iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)
    else
       ix = 0
       iy = 0
       scale = 1
       ismax = 4
       do i=0, ismax
          ip_low = iand(ipf,1023_MKD)
          ix = ix + scale * pix2x(ip_low)
          iy = iy + scale * pix2y(ip_low)
          scale = scale * 32
          ipf   = ipf/1024
       enddo
       ix = ix + scale * pix2x(ipf)
       iy = iy + scale * pix2y(ipf)
    endif

    !     transforms this in (horizontal, vertical) coordinates
    jrt = ix + iy  ! 'vertical' in {0,2*(nside-1)}
    jpt = ix - iy  ! 'horizontal' in {-nside+1,nside-1}

    !     computes the z coordinate on the sphere
    jr =  jrll(face_num+1)*nside - jrt - 1   ! ring number in {1,4*nside-1}

    if (jr < nside) then     ! north pole region
       nr = jr
       !z = 1.0_dp - nr * fact1 * nr
       theta = 2.0_dp * asin( nr / (sqrt(6.0_dp) * fn) )
       !kshift = 0

    else if (jr <= 3*nside) then ! equatorial region
       !fact2 = 2.0_dp/(3.0_dp*fn)
       nr = nside
       theta = ACOS((2*nside-jr)* 2.0_dp/(3.0_dp*fn) )
       !kshift = iand(jr - nside, 1)

    else if (jr > 3*nside) then ! south pole region
       nr = nl4 - jr
       !z = - 1.0_dp + nr * fact1 * nr
       theta = PI - 2.0_dp * asin( nr / (sqrt(6.0_dp) * fn) )
       !kshift = 0
    endif


!     !     computes the phi coordinate on the sphere, in [0,2Pi]
!     jp = (jpll(face_num+1)*nr + jpt + 1_MKD + kshift)/2  ! 'phi' number in the ring in {1,4*nr}
!     if (jp > nl4) jp = jp - nl4
!     if (jp < 1)   jp = jp + nl4

!     phi = (jp - (kshift+1)*0.5_dp) * (halfpi / nr)
    !     computes the phi coordinate on the sphere, in [0,2Pi]
    jp = jpll(face_num+1)*nr + jpt  ! 'phi' number in the ring in {0,8*nr-1}
    if (jp < 0)   jp = jp + 2_MKD*nl4

    phi = jp  * (quartpi / nr)

    return

#ifdef DOI8B
  end subroutine pix2ang_nest_8
#else
  end subroutine pix2ang_nest
#endif
!=======================================================================
!     pix2vec_nest
!
!     renders vector (x,y,z) coordinates of the nominal pixel center
!     for the pixel number ipix (NESTED scheme)
!     given the map resolution parameter nside
!     also returns the (x,y,z) position of the 4 pixel vertices (=corners)
!     in the order N,W,S,E
!=======================================================================
#ifdef DOI8B
  subroutine pix2vec_nest_8(nside, ipix, vector, vertex)
    integer(i4b), parameter :: MKD = i8b
#else
  subroutine pix2vec_nest  (nside, ipix, vector, vertex)
    integer(i4b), parameter :: MKD = i4b
#endif
    INTEGER(KIND=I4B), INTENT(IN) :: nside
    INTEGER(KIND=MKD), INTENT(IN) :: ipix
    REAL(KIND=DP),     INTENT(OUT), dimension(1:) :: vector
    REAL(KIND=DP),     INTENT(OUT), dimension(1:,1:), optional :: vertex

    INTEGER(KIND=MKD) :: npix, npface, ipf
    INTEGER(KIND=I4B) :: ip_low, ip_trunc, ip_med, ip_hi
    INTEGER(KIND=I4B) :: face_num, ix, iy, kshift, scale, i, ismax
    INTEGER(KIND=I4B) :: jrt, jr, nr, jpt, jp, nl4
    REAL(KIND=DP)     :: z, fn, fact1, fact2, sth, phi

    ! coordinate of the lowest corner of each face
    INTEGER(KIND=I4B), dimension(1:12) :: jrll = (/ 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 /) ! in unit of nside
    INTEGER(KIND=I4B), dimension(1:12) :: jpll = (/ 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 /) ! in unit of nside/2

    real(kind=DP) :: phi_nv, phi_wv, phi_sv, phi_ev, phi_up, phi_dn, sin_phi, cos_phi
    real(kind=DP) :: z_nv, z_sv, sth_nv, sth_sv
    real(kind=DP) :: hdelta_phi
    integer(kind=I4B) :: iphi_mod, iphi_rat
    logical(kind=LGT) :: do_vertex
    integer(kind=i4b) :: diff_phi
    character(len=*), parameter :: code = "pix2vec_nest"

    !-----------------------------------------------------------------------
#ifdef DOI8B
    if (nside <= ns_max4) then ! use faster 32-bit routine whenever possible
       if (present(vertex)) then
          call pix2vec_nest(nside, int(ipix, kind=i4b), vector, vertex)
       else
          call pix2vec_nest(nside, int(ipix, kind=i4b), vector)
       endif
       return
    endif
#else
    if (nside > ns_max4) call fatal_error(code//"> nside out of range")
#endif
    npix = nside2npix(nside)       ! total number of points
    if (ipix <0 .or. ipix>npix-1) call fatal_error(code//"> ipix out of range")

    !     initiates the array for the pixel number -> (x,y) mapping
    if (pix2x(1023) <= 0) call mk_pix2xy()

    npface = nside * int(nside, kind=MKD)
    nl4    = 4*nside

    !     finds the face, and the number in the face
    face_num = ipix/npface  ! face number in {0,11}
    ipf = MODULO(ipix,npface)  ! pixel number in the face {0,npface-1}

    do_vertex = .false.
    if (present(vertex)) then
       if (size(vertex,dim=1) >= 3 .and. size(vertex,dim=2) >= 4) then
          do_vertex = .true.
       else
          call fatal_error(code//">  vertex array has wrong size ")
       endif
    endif
    fn = real(nside, kind=dp)
    fact1 = 1.0_dp/(3.0_dp*fn*fn)
    fact2 = 2.0_dp/(3.0_dp*fn)

    !     finds the x,y on the face (starting from the lowest corner)
    !     from the pixel number
    if (nside <= ns_max4) then
       ip_low = iand(ipf,1023_MKD)       ! content of the last 10 bits
       ip_trunc =    ipf/1024        ! truncation of the last 10 bits
       ip_med = iand(ip_trunc,1023)  ! content of the next 10 bits
       ip_hi  =      ip_trunc/1024   ! content of the high weight 10 bits

       ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
       iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)
    else
       ix = 0
       iy = 0
       scale = 1
       ismax = 4
       do i=0, ismax
          ip_low = iand(ipf,1023_MKD)
          ix = ix + scale * pix2x(ip_low)
          iy = iy + scale * pix2y(ip_low)
          scale = scale * 32
          ipf   = ipf/1024
       enddo
       ix = ix + scale * pix2x(ipf)
       iy = iy + scale * pix2y(ipf)
    endif

    !     transforms this in (horizontal, vertical) coordinates
    jrt = ix + iy  ! 'vertical' in {0,2*(nside-1)}
    jpt = ix - iy  ! 'horizontal' in {-nside+1,nside-1}

    !     computes the z coordinate on the sphere
    jr =  jrll(face_num+1)*nside - jrt - 1   ! ring number in {1,4*nside-1}


    if (jr < nside) then     ! north pole region
       nr = jr
       z = 1.0_dp - nr*fact1*nr
       sth = nr * sqrt(fact1 * (1.d0 + z) ) ! more accurate close to pole
       kshift = 0
       if (do_vertex) then
          z_nv = 1.0_dp - (nr-1)*fact1*(nr-1)
          z_sv = 1.0_dp - (nr+1)*fact1*(nr+1)
       endif

    else if (jr <= 3*nside) then ! equatorial region
       nr = nside
       z  = (2*nside-jr)*fact2
       sth = SQRT((1.0_dp-z)*(1.0_dp+z)) ! good enough on Equator
       kshift = iand(jr - nside, 1)
       if (do_vertex) then
          z_nv = (2*nside-jr+1)*fact2
          z_sv = (2*nside-jr-1)*fact2
          if (jr == nside) then ! northern transition
             z_nv =  1.0_dp - (nside-1) * fact1 * (nside-1)
          elseif (jr == 3*nside) then  ! southern transition
             z_sv = -1.0_dp + (nside-1) * fact1 * (nside-1)
          endif
       endif

    else if (jr > 3*nside) then ! south pole region
       nr = nl4 - jr
       z   = - 1.0_dp + nr*fact1*nr
       sth = nr * sqrt(fact1 * (1.d0 - z) )
       kshift = 0
       if (do_vertex) then
          z_nv = - 1.0_dp + (nr+1)*fact1*(nr+1)
          z_sv = - 1.0_dp + (nr-1)*fact1*(nr-1)
       endif
    endif

    !     computes the phi coordinate on the sphere, in [0,2Pi]
    jp = (jpll(face_num+1)*nr + jpt + 1_MKD + kshift)/2  ! 'phi' number in the ring in {1,4*nr}
    if (jp > nl4) jp = jp - nl4
    if (jp < 1)   jp = jp + nl4

    phi = (jp - (kshift+1)*0.5_dp) * (halfpi / nr)

    ! pixel center
    !
    cos_phi = cos(phi)
    sin_phi = sin(phi)
    vector(1) = sth * cos_phi
    vector(2) = sth * sin_phi
    vector(3) = z

    if (do_vertex) then
       phi_nv = phi
       phi_sv = phi
       diff_phi = 0 ! phi_nv = phi_sv = phi

       phi_up = 0.0_dp
       iphi_mod = MODULO(jp-1, nr) ! in {0,1,... nr-1}
       iphi_rat = (jp-1) / nr      ! in {0,1,2,3}
       if (nr > 1) phi_up = HALFPI * (iphi_rat +  iphi_mod   /real(nr-1,kind=dp))
       phi_dn             = HALFPI * (iphi_rat + (iphi_mod+1)/real(nr+1,kind=dp))
       if (jr < nside) then            ! North polar cap
          phi_nv = phi_up
          phi_sv = phi_dn
          diff_phi = 3 ! both phi_nv and phi_sv different from phi
       else if (jr > 3*nside) then     ! South polar cap
          phi_nv = phi_dn
          phi_sv = phi_up
          diff_phi = 3 ! both phi_nv and phi_sv different from phi
       else if (jr == nside) then      ! North transition
          phi_nv = phi_up
          diff_phi = 1
       else if (jr == 3*nside) then    ! South transition
          phi_sv = phi_up
          diff_phi = 2
       endif

       hdelta_phi = PI / (4.0_dp*nr)

       ! west vertex
       phi_wv      = phi - hdelta_phi
       vertex(1,2) = sth * COS(phi_wv)
       vertex(2,2) = sth * SIN(phi_wv)
       vertex(3,2) = z

       ! east vertex
       phi_ev      = phi + hdelta_phi
       vertex(1,4) = sth * COS(phi_ev)
       vertex(2,4) = sth * SIN(phi_ev)
       vertex(3,4) = z

       ! north and south vertices
       sth_nv = SQRT((1.0_dp-z_nv)*(1.0_dp+z_nv))
       sth_sv = SQRT((1.0_dp-z_sv)*(1.0_dp+z_sv))
       if (diff_phi == 0) then
          vertex(1,1) = sth_nv * cos_phi
          vertex(2,1) = sth_nv * sin_phi
          vertex(1,3) = sth_sv * cos_phi
          vertex(2,3) = sth_sv * sin_phi
       else
          vertex(1,1) = sth_nv * COS(phi_nv)
          vertex(2,1) = sth_nv * SIN(phi_nv)
          vertex(1,3) = sth_sv * COS(phi_sv)
          vertex(2,3) = sth_sv * SIN(phi_sv)
       endif
       vertex(3,1) = z_nv
       vertex(3,3) = z_sv
    endif

    return

#ifdef DOI8B
  end subroutine pix2vec_nest_8
#else
  end subroutine pix2vec_nest
#endif

!=======================================================================
!     ang2pix_nest
!
!     renders the pixel number ipix (NESTED scheme) for a pixel which contains
!     a point on a sphere at coordinates theta and phi, given the map
!     resolution parameter nside
!
! 2009-03-09: calculations done directly at nside rather than ns_max
!             with ifort, for x and y integers>0, 
!             iand(x,y-1) is slightly faster than modulo(x,y)
!=======================================================================
#ifdef DOI8B
  subroutine ang2pix_nest_8(nside, theta, phi, ipix)
    integer(i4b), parameter :: MKD = i8b
#else
  subroutine ang2pix_nest  (nside, theta, phi, ipix)
    integer(i4b), parameter :: MKD = i4b
#endif
    INTEGER(KIND=I4B), INTENT(IN)  :: nside
    REAL(KIND=DP),     INTENT(IN)  :: theta, phi
    INTEGER(KIND=MKD), INTENT(OUT) :: ipix

    integer(kind=MKD) :: ipf, scale, scale_factor
    REAL(KIND=DP)     ::  z, za, tt, tp, tmp
    INTEGER(KIND=I4B) :: jp, jm, ifp, ifm, face_num, &
         &     ix, iy, ix_low, iy_low, ntt, i, ismax, ipix4
    character(len=*), parameter :: code = "ang2pix_nest"
    !-----------------------------------------------------------------------
#ifdef DOI8B
    if (nside <= ns_max4) then ! use faster 32-bit routine whenever possible
       call ang2pix_nest(nside, theta, phi, ipix4)
       ipix = int(ipix4, kind=MKD)
       return
    endif
    if (nside <1 .or. nside > ns_max ) call fatal_error(code//"> nside out of range")
#else
    if (nside <1 .or. nside > ns_max4) call fatal_error(code//"> nside out of range")
#endif
    if (theta<0.0_dp .or. theta>pi)  then
       print*,code//"> theta : ",theta," is out of range [0,Pi]"
       call fatal_error
    endif
    if (x2pix1(127) <= 0) call mk_xy2pix1()

    z  = COS(theta)
    za = ABS(z)
    tt = MODULO(phi, twopi) / halfpi  ! in [0,4[

    if (za <= twothird) then ! equatorial region

       !        (the index of edge lines increase when the longitude=phi goes up)
       jp = INT(nside*(0.5_dp + tt - z*0.75_dp)) !  ascending edge line index
       jm = INT(nside*(0.5_dp + tt + z*0.75_dp)) ! descending edge line index

       !        finds the face
       ifp = jp / nside  ! in {0,4}
       ifm = jm / nside
       if (ifp == ifm) then          ! faces 4 to 7
          face_num = iand(ifp,3) + 4
       else if (ifp < ifm) then     ! (half-)faces 0 to 3
          face_num = iand(ifp,3)
       else                            ! (half-)faces 8 to 11
          face_num = iand(ifm,3) + 8
       endif

       ix =         iand(jm, nside-1)
       iy = nside - iand(jp, nside-1) - 1

    else ! polar region, za > 2/3

       ntt = INT(tt)
       if (ntt >= 4) ntt = 3
       tp = tt - ntt
!        tmp = SQRT( 3.0_dp*(1.0_dp - za) )  ! in ]0,1]
       if (z > 0.0_dp) then
          tmp = sqrt(6.0_dp) * sin( theta * 0.5_dp)
       else
          tmp = sqrt(6.0_dp) * cos( theta * 0.5_dp)
       endif

       !        (the index of edge lines increase when distance from the closest pole goes up)
       jp = INT( nside * tp          * tmp ) ! line going toward the pole as phi increases
       jm = INT( nside * (1.0_dp - tp) * tmp ) ! that one goes away of the closest pole
       jp = MIN(nside-1, jp) ! for points too close to the boundary
       jm = MIN(nside-1, jm)

       !        finds the face and pixel's (x,y)
       if (z >= 0) then
          face_num = ntt  ! in {0,3}
          ix = nside - jm - 1
          iy = nside - jp - 1
       else
          face_num = ntt + 8 ! in {8,11}
          ix =  jp
          iy =  jm
       endif

    endif

    if (nside <= ns_max4) then 
       ix_low = iand(ix, 127)
       iy_low = iand(iy, 127)
       ipf =     x2pix1(ix_low) + y2pix1(iy_low) &
            & + (x2pix1(ix/128) + y2pix1(iy/128)) * 16384
    else
       scale = 1_MKD
       scale_factor = 16384_MKD ! 128*128
       ipf = 0_MKD
       ismax = 1 ! for nside in [2^14, 2^20]
       if (nside >  1048576 ) ismax = 3
       do i=0, ismax
          ix_low = iand(ix, 127) ! last 7 bits
          iy_low = iand(iy, 127) ! last 7 bits
          ipf = ipf + (x2pix1(ix_low)+y2pix1(iy_low)) * scale
          scale = scale * scale_factor
          ix  =     ix / 128 ! truncate out last 7 bits
          iy  =     iy / 128
       enddo
       ipf =  ipf + (x2pix1(ix)+y2pix1(iy)) * scale
    endif
    ipix = ipf + face_num* int(nside,MKD) * nside    ! in {0, 12*nside**2 - 1}

    return

#ifdef DOI8B
  end subroutine ang2pix_nest_8
#else
  end subroutine ang2pix_nest
#endif

!=======================================================================
!     vec2pix_nest
!
!     renders the pixel number ipix (NESTED scheme) for a pixel which contains
!     a point on a sphere at coordinate vector (=x,y,z), given the map
!     resolution parameter nside
!
! 2009-03-10: calculations done directly at nside rather than ns_max
!=======================================================================
#ifdef DOI8B
  subroutine vec2pix_nest_8(nside, vector, ipix)
    integer(i4b), parameter :: MKD = I8B
    integer(i4b), parameter :: mykind = I4B ! replace with i8b to reach nside=2^29
#else
  subroutine vec2pix_nest  (nside, vector, ipix)
    integer(i4b), parameter :: MKD = I4B
    integer(i4b), parameter :: mykind = I4B
#endif
    INTEGER(KIND=I4B), INTENT(IN)                :: nside
    REAL(KIND=DP),     INTENT(IN), dimension(1:) :: vector
    INTEGER(KIND=MKD), INTENT(OUT)               :: ipix

    integer(kind=MKD) :: ipf, scale, scale_factor
    REAL(KIND=DP)     ::  z, za, tt, tp, tmp, dnorm, phi
    INTEGER(KIND=I4B) ::  jp, jm, ifp, ifm, face_num, &
         &     ix, iy, ix_low, iy_low, ntt, i, ismax, ipix4
    character(len=*), parameter :: code = "vec2pix_nest"

    !-----------------------------------------------------------------------
#ifdef DOI8B
    if (nside <= ns_max4) then ! use faster 32-bit routine whenever possible
       call vec2pix_nest(nside, vector, ipix4)
       ipix = int(ipix4, kind=MKD)
       return
    endif
    if (nside<1 .or. nside>ns_max ) call fatal_error(code//"> nside out of range")
#else
    if (nside<1 .or. nside>ns_max4) call fatal_error(code//"> nside out of range")
#endif
    dnorm = SQRT(vector(1)**2+vector(2)**2+vector(3)**2)
    z = vector(3) / dnorm
    phi = 0.0_dp
    if (vector(1) /= 0.0_dp .or. vector(2) /= 0.0_dp) &
         &     phi = ATAN2(vector(2),vector(1)) ! phi in ]-pi,pi]

    za = ABS(z)
    if (phi < 0.0)    phi = phi + twopi ! phi in [0,2pi[
    tt = phi / halfpi ! in [0,4[
    if (x2pix1(127) <= 0) call mk_xy2pix1()

    if (za <= twothird) then ! equatorial region

       !        (the index of edge lines increase when the longitude=phi goes up)
       jp = INT(nside*(0.5_dp + tt - z*0.75_dp)) !  ascending edge line index
       jm = INT(nside*(0.5_dp + tt + z*0.75_dp)) ! descending edge line index

       !        finds the face
       ifp = jp / nside  ! in {0,4}
       ifm = jm / nside
       if (ifp == ifm) then          ! faces 4 to 7
          face_num = iand(ifp,3) + 4
       else if (ifp < ifm) then     ! (half-)faces 0 to 3
          face_num = iand(ifp,3)
       else                            ! (half-)faces 8 to 11
          face_num = iand(ifm,3) + 8
       endif

       ix =         iand(jm, nside-1)
       iy = nside - iand(jp, nside-1) - 1

    else ! polar region, za > 2/3

       ntt = INT(tt)
       if (ntt >= 4) ntt = 3
       tp = tt - ntt
       !tmp = SQRT( 3.0_dp*(1.0_dp - za) )  ! in ]0,1]
       tmp = SQRT(vector(1)**2+vector(2)**2) / dnorm ! sin(theta)
       tmp = tmp * SQRT( 3.0_dp / (1.0_dp + za) ) !more accurate

       !        (the index of edge lines increase when distance from the closest pole goes up)
       jp = INT( nside * tp          * tmp ) ! line going toward the pole as phi increases
       jm = INT( nside * (1.0_dp - tp) * tmp ) ! that one goes away of the closest pole
       jp = MIN(nside-1, jp) ! for points too close to the boundary
       jm = MIN(nside-1, jm)

       !        finds the face and pixel's (x,y)
       if (z >= 0) then
          face_num = ntt  ! in {0,3}
          ix = nside - jm - 1
          iy = nside - jp - 1
       else
          face_num = ntt + 8 ! in {8,11}
          ix =  jp
          iy =  jm
       endif

    endif

    if (nside <= ns_max4) then 
       ix_low = iand(ix, 127)
       iy_low = iand(iy, 127)
       ipf =     x2pix1(ix_low) + y2pix1(iy_low) &
            & + (x2pix1(ix/128) + y2pix1(iy/128)) * 16384
    else
       scale = 1_MKD
       scale_factor = 16384_MKD ! 128*128
       ipf = 0_MKD
       ismax = 1 ! for nside in [2^14, 2^20]
       if (nside >  1048576 ) ismax = 3
       do i=0, ismax
          ix_low = iand(ix, 127) ! last 7 bits
          iy_low = iand(iy, 127) ! last 7 bits
          ipf = ipf + (x2pix1(ix_low)+y2pix1(iy_low)) * scale
          scale = scale * scale_factor
          ix  =     ix / 128 ! truncate out last 7 bits
          iy  =     iy / 128
       enddo
       ipf =  ipf + (x2pix1(ix)+y2pix1(iy)) * scale
    endif
    ipix = ipf + face_num* int(nside,MKD) * nside    ! in {0, 12*nside**2 - 1}

    return
#ifdef DOI8B
  end subroutine vec2pix_nest_8
#else
  end subroutine vec2pix_nest
#endif

!=======================================================================
!     nest2ring
!
!     performs conversion from NESTED to RING pixel number
!=======================================================================
#ifdef DOI8B
  subroutine nest2ring_8(nside, ipnest, ipring)
    integer(i4b), parameter :: MKD = I8B
#else
  subroutine nest2ring  (nside, ipnest, ipring)
    integer(i4b), parameter :: MKD = I4B
#endif
    INTEGER(KIND=I4B), INTENT(IN)  :: nside
    INTEGER(KIND=MKD), INTENT(IN)  :: ipnest
    INTEGER(KIND=MKD), INTENT(OUT) :: ipring

    INTEGER(KIND=MKD) :: npix, npface, n_before, ipf
    INTEGER(KIND=I4B) :: ip_low, ip_trunc, ip_med, ip_hi
    INTEGER(KIND=I4B) :: face_num, ix, iy, kshift, scale, i, ismax
    INTEGER(KIND=I4B) :: jrt, jr, nr, jpt, jp, nl4
    integer(kind=i4b) :: ipring4, ipnest4
    character(len=*), parameter :: code = "nest2ring"

    !-----------------------------------------------------------------------
#ifdef DOI8B
    if (nside <= ns_max4) then ! use faster 32-bit routine whenever possible
       ipnest4 = ipnest
       call nest2ring(nside, ipnest4, ipring4)
       ipring = ipring4
       return
    endif
#else
    if (nside>ns_max4) call fatal_error(code//"> nside out of range")
#endif
    npix = nside2npix(nside)
    if (ipnest<0 .or. ipnest>npix-1) call fatal_error(code//"> ipnest out of range")

    !     initiates the array for the pixel number -> (x,y) mapping
    if (pix2x(1023) <= 0) call mk_pix2xy()

    npface = nside*int(nside, MKD)
    nl4    = 4*nside

    !     finds the face, and the number in the face
    face_num = ipnest/npface  ! face number in {0,11}
    ipf = iand(ipnest, npface-1)  ! pixel number in the face {0,npface-1}

    !     finds the x,y on the face (starting from the lowest corner)
    !     from the pixel number
    if (nside <= ns_max4) then
       ip_low = iand(ipf,1023_MKD)       ! content of the last 10 bits
       ip_trunc =    ipf/1024        ! truncation of the last 10 bits
       ip_med = iand(ip_trunc,1023)  ! content of the next 10 bits
       ip_hi  =      ip_trunc/1024   ! content of the high weight 10 bits

       ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
       iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)
    else
       ix = 0
       iy = 0
       scale = 1
       ismax = 4
       do i=0, ismax
          ip_low = iand(ipf,1023_MKD)
          ix = ix + scale * pix2x(ip_low)
          iy = iy + scale * pix2y(ip_low)
          scale = scale * 32
          ipf   = ipf/1024
       enddo
       ix = ix + scale * pix2x(ipf)
       iy = iy + scale * pix2y(ipf)
    endif

    !     transforms this in (horizontal, vertical) coordinates
    jrt = ix + iy  ! 'vertical' in {0,2*(nside-1)}
    jpt = ix - iy  ! 'horizontal' in {-nside+1,nside-1}

    !     computes the z coordinate on the sphere
    jr =  jrll1(face_num)*nside - jrt - 1   ! ring number in {1,4*nside-1}

    if (jr < nside) then     ! north pole region
       nr = jr
       n_before = 2 * nr * (nr - 1_MKD)
       kshift = 0

    else if (jr <= 3*nside) then  ! equatorial region (the most frequent)
       nr = nside                
       n_before = 2 * nr * ( 2 * jr - nr - 1_MKD)
       kshift = iand(jr - nside, 1)

    else ! south pole region
       nr = nl4 - jr
       n_before = npix - 2 * nr * (nr + 1_MKD)
       kshift = 0

    endif

    !     computes the phi coordinate on the sphere, in [0,2Pi]
    jp = (jpll1(face_num)*nr + jpt + 1_MKD + kshift)/2  ! 'phi' number in the ring in {1,4*nr}
    if (jp > nl4) jp = jp - nl4
    if (jp < 1)   jp = jp + nl4

    ipring = n_before + jp - 1_MKD ! in {0, npix-1}

    return

#ifdef DOI8B
  end subroutine nest2ring_8
#else
  end subroutine nest2ring
#endif

!=======================================================================
!     ring2nest
!
!     performs conversion from RING to NESTED pixel number
!=======================================================================
#ifdef DOI8B
  subroutine ring2nest_8(nside, ipring, ipnest)
    integer(i4b), parameter :: MKD = I8B
#else
  subroutine ring2nest  (nside, ipring, ipnest)
    integer(i4b), parameter :: MKD = I4B
#endif
    INTEGER(KIND=I4B), INTENT(IN)  :: nside
    INTEGER(KIND=MKD), INTENT(IN)  :: ipring
    INTEGER(KIND=MKD), INTENT(OUT) :: ipnest

    INTEGER(KIND=MKD) :: npix, ncap, ip
    INTEGER(KIND=MKD) :: scale, scale_factor
    INTEGER(KIND=I4B) :: nl2, nl4, iphi, ipt, &
         &     kshift, face_num, nr, &
         &     irn, ire, irm, irs, irt, ifm , ifp, &
         &     ix, iy, ix_low, iy_low, ipf, i, ismax
    integer(kind=i4b) :: ipring4, ipnest4
    character(len=*), parameter :: code = "ring2nest"

    !-----------------------------------------------------------------------
#ifdef DOI8B
    if (nside <= ns_max4) then ! use faster 32-bit routine whenever possible
       ipring4 = ipring
       call ring2nest(nside, ipring4, ipnest4)
       ipnest = ipnest4
       return
    endif
#else
    if (nside>ns_max4) call fatal_error(code//"> nside out of range")
#endif
    npix = nside2npix(nside)
    if (ipring <0 .or. ipring>npix-1) call fatal_error(code//"> ipring out of range")

    if (x2pix1(127) <= 0) call mk_xy2pix1()

    nl2 = 2*nside
    nl4 = 4*nside
    ncap = nl2*(nside-1_MKD) ! points in each polar cap, =0 for nside =1

    !     finds the ring number, the position of the ring and the face number
    if (ipring < ncap) then ! north polar cap

!       irn   = nint( sqrt( (ipring+1) * 0.50000_dp)) ! counted from North pole
       irn   = (cheap_isqrt(2*ipring+2) + 1) / 2 ! counted from North pole
       iphi  = ipring - 2*irn*(irn - 1_MKD)
! #ifdef DOI8B
!        ! fix round-off error appearing at large Nside
!        call correct_ring_phi(1, irn, iphi)
! #endif

       kshift = 0
       nr = irn                  ! 1/4 of the number of points on the current ring
       face_num = iphi / irn ! in {0,3}

    elseif (ipring < npix - ncap) then ! equatorial region

       ip    = ipring - ncap
       irn   = INT( ip / nl4 ) + nside               ! counted from North pole
       iphi  = iand(ip, nl4-1_MKD)

       kshift  = iand(irn+nside,1) ! MODULO(irn+nside,2)  ! 1 if irn+nside is odd, 0 otherwise
       nr = nside
       ire =  irn - nside + 1 ! in {1, 2*nside +1}
       irm =  nl2 + 2 - ire
       ifm = (iphi - ire/2 + nside) / nside ! face boundary
       ifp = (iphi - irm/2 + nside) / nside
       if (ifp == ifm) then          ! faces 4 to 7
          face_num = iand(ifp,3) + 4
       else if (ifp < ifm) then ! (half-)faces 0 to 3
          face_num = ifp
       else ! (half-)faces 8 to 11
          face_num = ifp + 7
       endif

    else ! south polar cap

       ip    = npix - ipring
!       irs   = nint( sqrt( ip * 0.500000_dp) )  ! counted from South pole
       irs   = (cheap_isqrt(2*ip) +1)/2  ! counted from South pole
       iphi  = 2*irs*(irs + 1_MKD) - ip
! #ifdef DOI8B
!        ! fix round-off error appearing at large Nside
!        call correct_ring_phi(-1, irs, iphi)
! #endif

       kshift = 0
       nr = irs
       irn   = nl4 - irs
       face_num = iphi / irs + 8 ! in {8,11}

    endif

    !     finds the (x,y) on the face
    irt =   irn  - jrll1(face_num)*nside + 1       ! in {-nside+1,0}
    ipt = 2*iphi - jpll1(face_num)*nr - kshift + 1 ! in {-nside+1,nside-1}
    if (ipt >= nl2) ipt = ipt - 8*nside ! for the face #4

    ix =  (ipt - irt ) / 2
    iy = -(ipt + irt ) / 2

    if (nside <= ns_max4) then 
       ix_low = iand(ix, 127)
       iy_low = iand(iy, 127)
       ipf =     x2pix1(ix_low) + y2pix1(iy_low) &
            & + (x2pix1(ix/128) + y2pix1(iy/128)) * 16384
    else
       scale = 1_MKD
       scale_factor = 16384_MKD ! 128*128
       ipf = 0_MKD
       ismax = 1 ! for nside in [2^14, 2^20]
       if (nside >  1048576 ) ismax = 3
       do i=0, ismax
          ix_low = iand(ix, 127) ! last 7 bits
          iy_low = iand(iy, 127) ! last 7 bits
          ipf = ipf + (x2pix1(ix_low)+y2pix1(iy_low)) * scale
          scale = scale * scale_factor
          ix  =     ix / 128 ! truncate out last 7 bits
          iy  =     iy / 128
       enddo
       ipf =  ipf + (x2pix1(ix)+y2pix1(iy)) * scale
    endif
    ipnest = ipf + face_num* (npix/12)    ! in {0, 12*nside**2 - 1}

    return

#ifdef DOI8B
  end subroutine ring2nest_8
#else
  end subroutine ring2nest
#endif

!=======================================================================
!   npix2nside
!
! given npix, returns nside such that npix = 12*nside^2
!  nside should be a power of 2 smaller than ns_max
!  if not, -1 is returned
! EH, Feb-2000
! 2009-03-05, edited, accepts 8-byte npix
!=======================================================================
#ifdef DOI8B
  function npix2nside_8(npix) result(nside_result)
    integer(i4b), parameter :: MKD = I8B
    INTEGER(KIND=MKD), PARAMETER  :: npix_max = (12_MKD*ns_max)*ns_max
#else
  function npix2nside  (npix) result(nside_result)
    integer(i4b), parameter :: MKD = I4B
    INTEGER(KIND=MKD), PARAMETER  :: npix_max = (12_MKD*ns_max4)*ns_max4
#endif
    INTEGER(KIND=MKD), INTENT(IN) :: npix
    INTEGER(KIND=MKD)             :: npix1, npix2
    INTEGER(KIND=I4B)             :: nside_result
    INTEGER(KIND=I4B)             :: nside
    CHARACTER(LEN=*),  PARAMETER  :: code = "npix2nside"
    !=======================================================================

    if (npix < 12 .or. npix > npix_max) then
       print*, code//"> Npix="//trim(string(npix))// &
            & " is out of allowed range: {12,"//trim(string(npix_max))//"}"
       nside = -1
       goto 1
    endif

    nside = NINT( sqrt(npix/12.0_dp) )
    npix1 = (12_MKD*nside)*nside
    if (abs(npix1-npix) > 0) then
       print*, code//"> Npix="//trim(string(npix))// &
            & " is not 12 * Nside * Nside "
       nside = -1
       goto 1
    endif

    ! test validity of Nside
    npix2 = nside2npix(nside)
    if (npix2 < 0) then
       nside = -1
       goto 1
    endif

1   continue
    nside_result = nside
    return

#ifdef DOI8B
  end function npix2nside_8
#else
  end function npix2nside
#endif


!=======================================================================
! in_ring
!     returns the list of pixels in RING or NESTED scheme (listir)
!     and their number (nir)
!     with latitude in [phi0-dphi, phi0+dphi] on the ring ir
!     (in {1,4*nside-1})
!     the pixel id-numbers are in {0,12*nside^2-1}
!     the indexing is RING, unless NEST is set to 1
!=======================================================================
#ifdef DOI8B
  subroutine in_ring_8(nside, iz, phi0, dphi, listir, nir, nest)
    integer(i4b), parameter :: MKD = I8B
#else
  subroutine in_ring (nside, iz, phi0, dphi, listir, nir, nest)
    integer(i4b), parameter :: MKD = I4B
#endif
    !=======================================================================
    integer(kind=i4b), intent(in)                 :: nside, iz
    integer(kind=i4b), intent(out)                :: nir
    real(kind=dp),     intent(in)                 :: phi0, dphi
    integer(kind=MKD), intent(out), dimension(0:) :: listir
    integer(kind=i4b), intent(in), optional       :: nest

!     logical(kind=lgt) :: conservative = .true.
    logical(kind=lgt) :: conservative = .false.
    logical(kind=lgt) :: take_all, to_top, do_ring

    integer(kind=i4b) :: i, diff
    integer(kind=i4b) :: nir1, nir2, ir, kshift
    integer(kind=MKD) :: ii, in, inext, npix, ncap, ipix1, ipix2
    integer(kind=MKD) :: ip_low, ip_hi, nr
    real(kind=dp)     :: phi_low, phi_hi, shift
    !=======================================================================

    take_all = .false.
    to_top   = .false.
    do_ring  = .true.
    if (present(nest)) then
       do_ring = (nest == 0)
    endif
    npix = (12_MKD * nside) * nside
    ncap  = 2*nside*(nside-1_MKD) ! number of pixels in the north polar cap
    listir = -1_MKD
    nir = 0

    phi_low = MODULO(phi0 - dphi, twopi)
    phi_hi  = MODULO(phi0 + dphi, twopi)
    if (ABS(dphi-PI) < 1.0e-6_dp) take_all = .true.

    !     ------------ identifies ring number --------------
    if (iz >= nside .and. iz <= 3*nside) then ! equatorial region
       ir = iz - nside + 1  ! in {1, 2*nside + 1}
       ipix1 = ncap  + 4*nside*(ir-1_MKD) !  lowest pixel number in the ring
       ipix2 = ipix1 + 4*nside - 1_MKD    ! highest pixel number in the ring
       kshift = MODULO(ir,2)
       nr = nside*4
    else
       if (iz < nside) then       !    north pole
          ir = iz
          ipix1 = 2*ir*(ir-1_MKD)        !  lowest pixel number in the ring
          ipix2 = ipix1 + 4*ir - 1_MKD   ! highest pixel number in the ring
       else                          !    south pole
          ir = 4*nside - iz
          ipix1 = npix - 2*ir*(ir+1_MKD) !  lowest pixel number in the ring
          ipix2 = ipix1 + 4*ir - 1_MKD   ! highest pixel number in the ring
       endif
       nr = ir*4
       kshift = 1
    endif

    !     ----------- constructs the pixel list --------------
    if (take_all) then
       nir    = ipix2 - ipix1 + 1
       if (do_ring) then
          listir(0:nir-1) = (/ (ii, ii=ipix1,ipix2) /)
       else
          call ring2nest(nside, ipix1, in)
          listir(0) = in
          do i=1,nir-1
             call next_in_line_nest(nside, in, inext)
             in = inext
             listir(i) = in
          enddo
       endif
       return
    endif

    shift = kshift * 0.5_dp
    if (conservative) then
       ! conservative : include every intersected pixels,
       ! even if pixel CENTER is not in the range [phi_low, phi_hi]
       ip_low = nint (nr * phi_low / TWOPI - shift, kind=MKD)
       ip_hi  = nint (nr * phi_hi  / TWOPI - shift, kind=MKD)
       ip_low = modulo (ip_low, nr) ! in {0,nr-1}
       ip_hi  = modulo (ip_hi , nr) ! in {0,nr-1}
    else
       ! strict : include only pixels whose CENTER is in [phi_low, phi_hi]
       ip_low = ceiling (nr * phi_low / TWOPI - shift, kind=MKD)
       ip_hi  = floor   (nr * phi_hi  / TWOPI - shift, kind=MKD)
!        if ((ip_low - ip_hi == 1) .and. (dphi*nr < PI)) then ! EH, 2004-06-01
       diff = modulo(ip_low - ip_hi, nr) ! in {-nr+1, nr-1} or {0,nr-1} ???
       if (diff < 0) diff = diff + nr    ! in {0,nr-1}
       if ((diff == 1) .and. (dphi*nr < PI)) then
          ! the interval is so small (and away from pixel center)
          ! that no pixel is included in it
          nir = 0
          return
       endif
!        ip_low = min(ip_low, nr-1) !  EH, 2004-05-28
!        ip_hi  = max(ip_hi , 0   )
       if (ip_low >= nr) ip_low = ip_low - nr
       if (ip_hi  <  0 ) ip_hi  = ip_hi  + nr
    endif
    !
    if (ip_low > ip_hi) to_top = .true.
    ip_low = ip_low + ipix1
    ip_hi  = ip_hi  + ipix1

    if (to_top) then
       nir1 = ipix2 - ip_low + 1
       nir2 = ip_hi - ipix1  + 1
       nir  = nir1 + nir2
       if (do_ring) then
          listir(0:nir1-1)   = (/ (ii, ii=ip_low, ipix2) /)
          listir(nir1:nir-1) = (/ (ii, ii=ipix1, ip_hi) /)
       else
          call ring2nest(nside, ip_low, in)
          listir(0) = in
          do i=1,nir-1
             call next_in_line_nest(nside, in, inext)
             in = inext
             listir(i) = in
          enddo
       endif
    else
       nir = ip_hi - ip_low + 1
       if (do_ring) then
          listir(0:nir-1) = (/ (ii, ii=ip_low, ip_hi) /)
       else
          call ring2nest(nside, ip_low, in)
          listir(0) = in
          do i=1,nir-1
             call next_in_line_nest(nside, in, inext)
             in = inext
             listir(i) = in
          enddo
       endif
    endif

    return
#ifdef DOI8B
  end subroutine in_ring_8
#else
  end subroutine in_ring
#endif

  !=======================================================================
#ifdef DOI8B
  subroutine getdisc_ring_8( nside, vector0, radius, listpix, nlist)
    integer(i4b), parameter :: MKD = I8B
#else
  subroutine getdisc_ring ( nside, vector0, radius, listpix, nlist)
    integer(i4b), parameter :: MKD = I4B
#endif

    integer(kind=I4B), intent(in)                 :: nside
    real(kind=DP),     intent(in), dimension(1:)  :: vector0
    real(kind=DP),     intent(in)                 :: radius
    integer(kind=MKD), intent(out), dimension(0:) :: listpix
    integer(kind=MKD), intent(out)                :: nlist
  !=======================================================================

    print*,"-------------------------------------------------------------"
    print*,"WARNING : the routine getdisc_ring is now obsolete"
    print*,"  Use "
    print*," call query_disc(nside, vector0, radius_radian, listpix, nlist [, nest]) "
    print*,"  instead (same module)"
    print*," It lets you choose the indexing scheme (RING or NESTED) "
    print*," used for the outuput"
    print*,"-------------------------------------------------------------"

    call query_disc(nside, vector0, radius, listpix, nlist)

    return

#ifdef DOI8B
  end subroutine getdisc_ring_8
#else
  end subroutine getdisc_ring
#endif

!=======================================================================
#ifdef DOI8B
  subroutine query_disc_old_8( nside, vector0, radius, listpix, nlist, nest, inclusive)
    integer(i4b), parameter :: MKD = I8B
#else
  subroutine query_disc_old ( nside, vector0, radius, listpix, nlist, nest, inclusive)
    integer(i4b), parameter :: MKD = I4B
#endif
    !=======================================================================

    integer(kind=I4B), intent(in)                 :: nside
    real(kind=DP),     intent(in), dimension(1:)  :: vector0
    real(kind=DP),     intent(in)                 :: radius
    integer(kind=MKD), intent(out), dimension(0:) :: listpix
    integer(kind=MKD), intent(out)                :: nlist
    integer(kind=I4B), intent(in), optional       :: nest
    integer(kind=I4B), intent(in), optional       :: inclusive

    INTEGER(kind=I4B) :: irmin, irmax, iz, ip, nir
    REAL(kind=DP) :: norm_vect0
    REAL(kind=DP) :: x0, y0, z0, radius_eff, fudge
    REAL(kind=DP) :: a, b, c, cosang
    REAL(kind=DP) :: dth1, dth2
    REAL(kind=DP) :: phi0, cosphi0, cosdphi, dphi
    REAL(kind=DP) :: rlat0, rlat1, rlat2, zmin, zmax, z
    INTEGER(kind=MKD), DIMENSION(:),   ALLOCATABLE  :: listir
    integer(kind=MKD) :: npix, list_size, nlost, ilist
    INTEGER(kind=I4B) :: status
    character(len=*), parameter :: code = "QUERY_DISC_OLD"
    logical(kind=LGT) :: do_inclusive
    integer(kind=I4B)                             :: my_nest

    !=======================================================================

    list_size = long_size(listpix)
    !     ---------- check inputs ----------------
    npix = nside2npix(nside)

    if (radius < 0.0_dp .or. radius > PI) then
       write(unit=*,fmt="(a)") code//"> the angular radius is in RADIAN "
       write(unit=*,fmt="(a)") code//"> and should lie in [0,Pi] "
       call fatal_error("> program abort ")
    endif

    do_inclusive = .false.
    if (present(inclusive)) then
       if (inclusive == 1) do_inclusive = .true.
    endif

    my_nest = 0
    if (present(nest)) then
       if (nest == 0 .or. nest == 1) then
          my_nest = nest
       else
          print*,code//"> NEST should be 0 or 1"
          call fatal_error("> program abort ")
       endif
    endif

    !     --------- allocate memory -------------
    ALLOCATE( listir(0: 4*nside-1), STAT = status)
    if (status /= 0) then
       write(unit=*,fmt="(a)") code//"> can not allocate memory for listir :"
       call fatal_error("> program abort ")
    endif

    dth1 = 1.0_dp / (3.0_dp*real(nside,kind=dp)**2)
    dth2 = 2.0_dp / (3.0_dp*real(nside,kind=dp))

    radius_eff = radius
    if (do_inclusive) then
! !        fudge = PI / (4.0_dp*nside) ! increase radius by half pixel size
!        fudge = acos(TWOTHIRD) / real(nside,kind=dp) ! 1.071* half pixel size
!        radius_eff = radius + fudge
       radius_eff = fudge_query_radius(nside, radius)
    endif
    cosang = COS(radius_eff)

    !     ---------- circle center -------------
    norm_vect0 =  SQRT(DOT_PRODUCT(vector0,vector0))
    x0 = vector0(1) / norm_vect0
    y0 = vector0(2) / norm_vect0
    z0 = vector0(3) / norm_vect0

    phi0=0.0_dp
    if ((x0/=0.0_dp).or.(y0/=0.0_dp)) phi0 = ATAN2 (y0, x0)  ! in ]-Pi, Pi]
    cosphi0 = COS(phi0)
    a = x0*x0 + y0*y0

    !     --- coordinate z of highest and lowest points in the disc ---
    rlat0  = ASIN(z0)    ! latitude in RAD of the center
    rlat1  = rlat0 + radius_eff
    rlat2  = rlat0 - radius_eff
    if (rlat1 >=  halfpi) then
       zmax =  1.0_dp
    else
       zmax = SIN(rlat1)
    endif
    irmin = ring_num(nside, zmax)
    irmin = MAX(1, irmin - 1) ! start from a higher point, to be safe

    if (rlat2 <= -halfpi) then
       zmin = -1.0_dp
    else
       zmin = SIN(rlat2)
    endif
    irmax = ring_num(nside, zmin)
    irmax = MIN(4*nside-1, irmax + 1) ! go down to a lower point

    ilist = -1_MKD

    !     ------------- loop on ring number ---------------------
    do iz = irmin, irmax

       z = ring2z(nside, iz) 

       !        --------- phi range in the disc for each z ---------
       b = cosang - z*z0
       c = 1.0_dp - z*z
       if ((x0==0.0_dp).and.(y0==0.0_dp)) then
          dphi=PI
          if (b > 0.0_dp) goto 1000 ! out of the disc, 2008-03-30
          goto 500
       endif
       cosdphi = b / SQRT(a*c)
       if (ABS(cosdphi) <= 1.0_dp) then
          dphi = ACOS (cosdphi) ! in [0,Pi]
       else
          if (cosphi0 < cosdphi) goto 1000 ! out of the disc
          dphi = PI ! all the pixels at this elevation are in the disc
       endif
500    continue

       !        ------- finds pixels in the disc ---------
       call in_ring(nside, iz, phi0, dphi, listir, nir, nest=my_nest)

       !        ----------- merge pixel lists -----------
       nlost = ilist + nir + 1_MKD - list_size
       if ( nlost > 0 ) then
          print*,code//"> listpix is too short, it will be truncated at ",nir
          print*,"                         pixels lost : ", nlost
          nir = nir - nlost
       endif
       do ip = 0, nir-1
          ilist = ilist + 1_MKD
          listpix(ilist) = listir(ip)
       enddo

1000   continue
    enddo

    !     ------ total number of pixel in the disc --------
    nlist = ilist + 1_MKD


    !     ------- deallocate memory and exit ------
    DEALLOCATE(listir)

    return

#ifdef DOI8B
  end subroutine query_disc_old_8
#else
  end subroutine query_disc_old
#endif


  !=======================================================================
#ifdef DOI8B
  subroutine discedge2fulldisc_8( nside, ringphi, ngr, list, nlist)
    integer(i4b), parameter :: MKD = I8B
#else
  subroutine discedge2fulldisc( nside, ringphi, ngr, list, nlist)
    integer(i4b), parameter :: MKD = I4B
#endif
    !=======================================================================
    integer(i4b), intent(in) :: nside
    integer(i4b), dimension(1:3,1:ngr), intent(in) :: ringphi
    integer(i4b),                       intent(in) :: ngr
    integer(MKD), dimension(1:),        intent(out) :: list
    integer(MKD),                       intent(out) :: nlist
    !
    integer(i4b) :: j, jr_min, jr_max, nj, nr, kshift, np, my_low, my_hi
    integer(i4b) :: ir, i, ip
    integer(i8b) :: npc, listsize
    !=======================================================================
    
    listsize = long_size(list)
    nlist = 0_MKD
    if (ngr == 0) then ! no valid rings
       list(1) = -1
       return
    endif

    do j=0, ngr-1
       ir = ringphi(1, j+1)             ! current ring, in [1, nl4-1]
       call pixels_per_ring(nside, ir, nr, kshift, npc)
       my_low = ringphi(2, j+1)
       if (my_low >= 0) then
          my_hi  = ringphi(3, j+1)
          np = my_hi - my_low     ! in [-nr+1, nr-1]
          np = modulo(np, nr) + 1 ! deal with periodic BC
          np = min(np, nr)
          if (nlist + np > listsize) then
             print*,'Pixel query: too many pixels found for output list provided.'
             print*,'truncated at ',nlist
             return
          endif
          do i=0, np-1
             ip = modulo(my_low + i, nr)
             list(nlist+1+i) = npc - nr + ip ! fill final list
          enddo
          nlist = nlist + np
       endif
    enddo
    
    if (nlist == 0) list(1) = -1
    
    return
#ifdef DOI8B
  end subroutine discedge2fulldisc_8
#else
  end subroutine discedge2fulldisc
#endif

!=======================================================================
!
!      query_disc (Nside, Vector0, Radius, Listpix, Nlist[, Nest, Inclusive])
!      ----------
!      routine for pixel query in the RING or NESTED scheme
!      all pixels within an angular distance Radius of the center
!
!     Nside    = resolution parameter (a power of 2)
!     Vector0  = central point vector position (x,y,z in double precision)
!     Radius   = angular radius in RADIAN (in double precision)
!     Listpix  = list of pixel closer to the center (angular distance) than Radius
!     Nlist    = number of pixels in the list
!     nest  (OPT), :0 by default, the output list is in RING scheme
!                  if set to 1, the output list is in NESTED scheme
!     inclusive (OPT) , :0 by default, only the pixels whose center
!                       lie in the triangle are listed on output
!                  if set to 1, all pixels overlapping the triangle are output
!
!      * all pixel numbers are in {0, 12*Nside*Nside - 1}
!     NB : the dimension of the listpix array is fixed in the calling
!     routine and should be large enough for the specific configuration
!
!      lower level subroutines called by getdisc_ring :
!       (you don't need to know them)
!      x=fudge_query_radius()
!      x=ring_num (nside, ir)
!      x=ring2z()
!      discphirange_at_z()
!      pixels_on_edge()
!      check_edge_pixels()
!      discedge2fulldisc()
!      -------
!
! v1.0, EH, TAC, ??
! v1.1, EH, Caltech, Dec-2001
! v1.2, EH, IAP, 2008-03-30: fixed bug appearing when disc centered on either pole
! 2009-06-17: deals with Nside > 8192
! 2011-06-09: uses ring2z
! 2011-10-18: improve fudge radius determination.
! New algorithm for Inclusive case: test boundary of edge pixels on each ring
!    2013-04-02: bug correction in query_disc in inclusive mode
!    2014-07-07: bug correction in discedge2fulldisc
! 2020-08-31: bug correction in check_edge_pixels (reported by R. de Belsunce)
    !=======================================================================
#ifdef DOI8B
  subroutine query_disc_8( nside, vector0, radius, listpix, nlist, nest, inclusive)
    integer(i4b), parameter :: MKD = I8B
#else
  subroutine query_disc ( nside, vector0, radius, listpix, nlist, nest, inclusive)
    integer(i4b), parameter :: MKD = I4B
#endif
    !=======================================================================
    integer(kind=I4B), intent(in)                 :: nside
    real(kind=DP),     intent(in), dimension(1:)  :: vector0
    real(kind=DP),     intent(in)                 :: radius
    integer(kind=MKD), intent(out), dimension(0:) :: listpix
    integer(kind=MKD), intent(out)                :: nlist
    integer(kind=I4B), intent(in), optional       :: nest
    integer(kind=I4B), intent(in), optional       :: inclusive

    INTEGER(kind=I4B) :: irmin, irmax, iz, ip, nir, nr, ngr, nrh
    REAL(kind=DP) :: norm_vect0
    REAL(kind=DP) :: z0, radius_eff, fudge
    REAL(kind=DP) :: phi0, dphi
    REAL(kind=DP) :: rlat0, rlat1, rlat2, zmin, zmax, z
    integer(kind=MKD) :: npix, list_size, nlost, ilist
    INTEGER(kind=I4B) :: status
    character(len=*), parameter :: code = "QUERY_DISC"
    logical(kind=LGT) :: do_inclusive
    integer(kind=I4B) :: my_nest
!    real(dp), allocatable, dimension(:) :: ztab, dphitab
    real(dp), dimension(1:4*nside-1) :: ztab, dphitab
    real(dp), dimension(:), allocatable :: zlist, dphilist
    integer(i4b), dimension(1:3, 1:4*nside-1) :: ringphi
    integer(i4b) :: nsideh, nsboost
    real(dp) :: radiush

    !=======================================================================

    list_size = long_size(listpix)
    !     ---------- check inputs ----------------
    npix = nside2npix(nside)

    if (radius < 0.0_dp .or. radius > PI) then
       write(unit=*,fmt="(a)") code//"> the angular radius is in RADIAN "
       write(unit=*,fmt="(a)") code//"> and should lie in [0,Pi] "
       call fatal_error("> program abort ")
    endif

    do_inclusive = .false.
    if (present(inclusive)) then
       if (inclusive == 1) do_inclusive = .true.
    endif

    my_nest = 0
    if (present(nest)) then
       if (nest == 0 .or. nest == 1) then
          my_nest = nest
       else
          print*,code//"> NEST should be 0 or 1"
          call fatal_error("> program abort ")
       endif
    endif

    radius_eff = radius
    if (do_inclusive) radius_eff = fudge_query_radius(nside, radius)

    !     ---------- circle center -------------
    norm_vect0 =  SQRT(DOT_PRODUCT(vector0,vector0))
    z0 = vector0(3) / norm_vect0

    !     --- coordinate z of highest and lowest points in the disc ---
    rlat0  = ASIN(z0)    ! latitude in RAD of the center
    rlat1  = rlat0 + radius_eff
    rlat2  = rlat0 - radius_eff
    if (rlat1 >=  halfpi) then
       zmax =  1.0_dp
    else
       zmax = SIN(rlat1)
    endif
    irmin = ring_num(nside, zmax)
    irmin = MAX(1, irmin - 1) ! start from a higher point, to be safe

    if (rlat2 <= -halfpi) then
       zmin = -1.0_dp
    else
       zmin = SIN(rlat2)
    endif
    irmax = ring_num(nside, zmin)
    irmax = MIN(4*nside-1, irmax + 1) ! go down to a lower point

    nr = irmax-irmin+1 ! in [1, 4*Nside-1]
    do iz = irmin, irmax
       ztab(iz-irmin+1) = ring2z(nside, iz)
    enddo
    call discphirange_at_z(vector0, radius_eff, ztab, nr, dphitab, phi0)
    call pixels_on_edge(nside, irmin, irmax, phi0, dphitab, ringphi, ngr)
    if (do_inclusive) then
       ! sample edge pixels at larger Nside
       nsboost = 2**ceiling( log(1.1_DP/(radius*nside))/log(2.0_DP) ) ! make pixel smaller than disc (2020-08-31)
       nsboost = max(nsboost,16)
       nsideh  = min(NS_MAX8, nside * int(nsboost,i8b))
       radiush = fudge_query_radius(nsideh, radius, quadratic=.true.)

       irmin = ring_num(nsideh, zmax, shift=+1) ! shifted South
       irmax = ring_num(nsideh, zmin, shift=-1) ! shifted North
       nrh = irmax - irmin + 1
       allocate(zlist(1:nrh), dphilist(1:nrh))
       do iz = irmin, irmax
          zlist(iz-irmin+1) = ring2z(nsideh, iz)
       enddo
       call discphirange_at_z(vector0, radiush, zlist, nrh, dphilist, phi0)
       phi0 = mod(phi0 + TWOPI, TWOPI) ! [-Pi,Pi] -> [0,2Pi]
       call check_edge_pixels(nside, nsboost, irmin, irmax, phi0, dphilist, ringphi, ngr)
       deallocate(zlist, dphilist)
    endif
    call discedge2fulldisc(nside, ringphi, ngr, listpix, nlist)

    if (my_nest == 1) then
       do ip=0_MKD, nlist-1
          call ring2nest(nside, listpix(ip), listpix(ip))
       enddo
    endif
    
    return

#ifdef DOI8B
  end subroutine query_disc_8
#else
  end subroutine query_disc
#endif

!=======================================================================
! query_strip ( nside, theta1, theta2, listpix, nlist, nest, inclusive)
!
! finds pixels having a colatitude (measured from North Pole):
!  theta1 < colatitude < theta2
!  with 0 <= theta1 < theta2 <= Pi
! if theta2 < theta1 then pixels with
! 0 <= colatitude < theta2   or   theta1 < colatitude < Pi are returned
!
! nside :          I4,       input,  resolution parameter
! theta1, theta2 : DP,       input,  bounds
! listpix :        I4/I8  array, output,  list of pixels
! nlist :          I4/I8,        output, number of pixels in list
! nest  :          I4 opt.,  input, =0 : RING scheme, =1 : NESTED scheme
! inclusive:       I4 opt.,  input, if=1 include pixel touched by strip, 
!                                   if 0 keep pixels with center in strip
!
! v1.0, EH, Caltech, Jan-2002
! v2.0: 2004-12 : added inclusive keyword
! v2.1: 2008-03-28: improved inclusive regime
! 2009-06-15: deals with Nside > 8192
!=======================================================================
#ifdef DOI8B
  subroutine query_strip_8( nside, theta1, theta2, listpix, nlist, nest, inclusive)
    integer(kind=i4b), parameter  ::   MKD = I8B
#else
  subroutine query_strip  ( nside, theta1, theta2, listpix, nlist, nest, inclusive)
    integer(kind=i4b), parameter  ::   MKD = I4B
#endif
    integer(kind=I4B), intent(in)                    :: nside
    real(kind=DP),     intent(in)                    :: theta1, theta2
    integer(kind=MKD), intent(out), dimension(0:)    :: listpix
    integer(kind=MKD), intent(out)                   :: nlist
    integer(kind=I4B), intent(in),  optional         :: nest
    integer(kind=i4b), intent(in),  optional         :: inclusive

    integer(kind=I4B)                  :: nstrip, my_nest
    integer(kind=I4B)                  :: iz, ip, is, irmin, irmax
    integer(kind=I4B)                  :: nir
    real(kind=DP)                      :: phi0, dphi, zu, zd, zring
    real(kind=DP), dimension(1:4)      :: colrange
    character(len=*), parameter        :: code = "query_strip"
    integer(kind=MKD)                  :: npix, ilist, nlost, list_size
    integer(kind=MKD), dimension(:), allocatable :: listir
    logical(kind=LGT) :: be_inclusive

    !=======================================================================
    list_size = long_size(listpix)
    !     ---------- check inputs ----------------
    npix = nside2npix(nside)
    if (npix < 0) then
       print*,code//"> Nside should be a power of 2"
       print*,code//"> current value = ",nside
       call fatal_error("> program abort ")
    endif

    if    (theta1 < 0.0_dp .or. theta1 > PI .or. &
         & theta2 < 0.0_dp .or. theta2 > PI) then
       write(unit=*,fmt="(a)") code//"> the colatitudes are in RADIAN "
       write(unit=*,fmt="(a)") code//"> and should lie in [0,Pi] "
       print*,code//"> current value = ", theta1, theta2
       call fatal_error("> program abort ")
    endif

    if (theta1 <= theta2) then
       nstrip = 1
       colrange(1:2*nstrip) = (/ theta1, theta2 /)
    else
       nstrip = 2
       colrange(1:2*nstrip) = (/ 0.0_dp, theta2, theta1, PI/)
    endif

    be_inclusive = .false.
    if (present(inclusive)) be_inclusive = (inclusive==1)

    my_nest = 0
    if (present(nest)) then
       if (nest == 0 .or. nest == 1) then
          my_nest = nest
       else
          print*,code//"> NEST should be 0 or 1"
          call fatal_error("> program abort ")
       endif
    endif

    !     ------------- loop on strips ---------------------
    ilist = -1_MKD
    allocate(listir(0:4*nside-1))
    do is =0, nstrip-1
       zu = cos(colrange(2*is+1)) ! upper bound in z
       zd = cos(colrange(2*is+2)) ! lower bound in z
       if (be_inclusive) then
          irmin = ring_num(nside, zu, shift=-1_i4b) ! shift up north
          irmax = ring_num(nside, zd, shift=+1_i4b) ! shift down south
       else
          irmin = ring_num(nside, zu)
          irmax = ring_num(nside, zd)
       endif


       !     ------------- loop on ring number ---------------------
       do iz = irmin, irmax

          zring = ring2z(nside, iz) ! z of ring being considered

          if ((zring >= zd .and. zring <= zu) .or. be_inclusive) then
             !        ------- finds pixels in the ring ---------
             phi0 = 0
             dphi = PI
             call in_ring(nside, iz, phi0, dphi, listir, nir, nest=my_nest)
             
             nlost = ilist + nir + 1_MKD - list_size
             if ( nlost > 0 ) then
                print*,code//"> listpix is too short, it will be truncated at ",nir
                print*,"                         pixels lost : ", nlost
                nir = nir - nlost
             endif
             do ip = 0, nir-1
                ilist = ilist + 1_MKD
                listpix(ilist) = listir(ip)
             enddo
          endif

       enddo

    enddo

    !     ------ total number of pixel in the strip --------
    nlist = ilist + 1_MKD

    !     ------- deallocate memory and exit ------
    DEALLOCATE(listir)

    return
#ifdef DOI8B
  end subroutine query_strip_8
#else
  end subroutine query_strip
#endif

!=======================================================================
!  query_polygon
!
! finds pixels that lie within a CONVEX polygon defined by its vertex on the sphere
!
! nside             : IN
! vlist(1:3, 0:n-1) : IN, list of vertices
! nv                : IN, number of vertices to be used (nv <= n)
! listpix           : OUT
! nlist             : OUT
! nest              : IN, OPTIONAL
! inclusive         : IN, OPTIONAL
!
! algorithm:
!   the polygon is divided into triangles
!   vertex 0 belongs to all the triangles
!
! v1.0, EH, Caltech, Dec-2001
! 2009-06-17: deals with Nside > 8192
  !=======================================================================
#ifdef DOI8B
  subroutine query_polygon_8( nside, vlist, nv, listpix, nlist, nest, inclusive)
    USE num_rec, ONLY : isort
    integer(kind=i4b), parameter  ::   MKD = I8B
#else
  subroutine query_polygon ( nside, vlist, nv, listpix, nlist, nest, inclusive)
    USE num_rec, ONLY : isort
    integer(kind=i4b), parameter  ::   MKD = I4B
#endif
    !=======================================================================
    integer(kind=I4B), intent(in)                    :: nside
    real(kind=DP),     intent(in),  dimension(1:,1:) :: vlist
    integer(kind=I4B), intent(in)                    :: nv
    integer(kind=MKD), intent(out), dimension(0:)    :: listpix
    integer(kind=MKD), intent(out)                   :: nlist
    integer(kind=I4B), intent(in),  optional         :: nest
    integer(kind=i4b), intent(in),  optional         :: inclusive

    real(kind=DP),     dimension(:), pointer         :: vp0, vp1, vp2
    real(kind=DP),     dimension(1:3)                :: vo
    real(kind=DP),     dimension(:,:), allocatable, target   :: vvlist
    real(kind=DP),     dimension(:),   allocatable   :: ss
    real(kind=DP)                                    :: surface, fsky, hand
    integer(kind=MKD)                                :: npix, n_in_trg, ilist, ntl
    integer(kind=MKD)                                :: list_size, nlost, i, k
    integer(kind=I4B)                                :: n_remain
    integer(kind=I4B)                                :: i0, i1, i2
    integer(kind=I4B)                                :: np, nm, nlow
    integer(kind=I4B), dimension(1:1)                :: ix
    integer(kind=MKD), dimension(:), allocatable     :: templist
    character(len=*), parameter :: code = "QUERY_POLYGON"
    integer(kind=I4B)                                :: my_nest

    !=======================================================================

    list_size = long_size(listpix)
    n_remain = nv
    ilist = -1_MKD
    listpix(:) = -1_MKD

    call assert(n_remain>=3, &
      code//"> The polygon should have at least 3 vertices")

    if (size(vlist) < n_remain*3) then
       print*,code//"> ",size(vlist)/3," vertices are given"
       print*,code//"> expected : ",n_remain
       call fatal_error
    endif

    allocate(vvlist(1:3, 1:n_remain))
    vvlist = vlist

    my_nest = 0
    if (present(nest)) then
       if (nest == 0 .or. nest == 1) then
          my_nest = nest
       else
          print*,code//"> NEST should be 0 or 1"
          call fatal_error("> program abort ")
       endif
    endif

    !-----------------------------------------------------------------
    ! check that the polygon is convex or has only one concave vertex
    !-----------------------------------------------------------------
    if (n_remain == 3) goto 1000 ! a triangle is always convex

    allocate(ss(1:n_remain))
    do i1 = 0, n_remain-1
       i0 = modulo(i1-1, n_remain)  ! in [0, n_remain-1]
       i2 = modulo(i1+1, n_remain)  ! in [0, n_remain-1]

       vp0 => vvlist(:,i0+1)
       vp1 => vvlist(:,i1+1)
       vp2 => vvlist(:,i2+1)

       ! computes the handedness   (v0 x v2) . v1  for each vertex v1
       call vect_prod(vp0, vp2, vo)
       hand = dot_product(vo, vp1)
       ss(i1+1) = sign(1.0_dp, hand) ! either +1 or -1
    enddo
    np = count( ss > 0.0_dp) ! number of vertex with positive handedness
    nm = n_remain - np

    nlow = min(np,nm)
    if (nlow == 0) goto 1000
    if (nlow == 1) then
       ! only one concave vertex
       if (np == 1) then
          ix = maxloc(ss)
       else
          ix = minloc(ss)
       endif
       ! rotate pixel list to put that vertex in #0
       vvlist = cshift(vvlist, ix(1)-1, dim = 2)
    endif
    if (nlow > 1) then
       ! more than 1 concave vertex
       print*,"***************************************"
       print*,"The polygon is not convex"
       print*," and has more than one concave vertex"
       print*,"The result is unpredictable"
       print*,"***************************************"
    endif
    deallocate(ss)

1000 continue

    !--------------------------------------------
    ! fill the polygon, one triangle at a time
    !--------------------------------------------
    npix = nside2npix(nside)
    do
       ! build triangle from vertices #0, n-2, and n-1
       vp0 => vvlist(1:3,1)
       vp1 => vvlist(1:3,n_remain-1)
       vp2 => vvlist(1:3,n_remain  )

       ! computes its surface
       call surface_triangle( vp0, vp1, vp2, surface )
       fsky = surface/FOURPI
       n_in_trg = npix * (fsky * 1.4) + 12*nside
       n_in_trg = min(n_in_trg, npix)

       ! find pixels within triangle
       allocate(templist(0:n_in_trg-1))
       call query_triangle( nside, vp0, vp1, vp2, templist, ntl, nest=my_nest,inclusive=inclusive )

       ! merge new list with existing one
       nlost = ilist + ntl + 1_MKD - list_size
       if ( nlost > 0 ) then
          print*,code//"> listpix is too short, it will be truncated at ",list_size
          print*,"                         pixels lost : ", nlost
          print*, list_size
          ntl = ntl - nlost
       endif
       do i = 0_MKD, ntl - 1_MKD
          ilist = ilist + 1_MKD
          listpix(ilist) = templist(i)
       enddo
       deallocate(templist)

       if (n_remain == 3) exit

       ! prune vertex list
       n_remain = n_remain - 1

    enddo
    deallocate(vvlist)

    !-------------------------
    ! make final pixel list
    !-------------------------
    nlist = ilist + 1_MKD

    ! sort final list
    if (nlist > 0) call isort(nlist, listpix)

    ! remove redondant pixels
    ! (we keep 0th element of the list)
    k = 1_MKD
    do i = 1, nlist-1_MKD
       if (listpix(i) > listpix(i-1)) then
          listpix(k) = listpix(i)
          k = k + 1_MKD
       endif
    enddo
    nlist = k

    ! pad end of list with -1
    listpix(nlist:list_size-1) = -1_MKD

    return
#ifdef DOI8B
  end subroutine query_polygon_8
#else
  end subroutine query_polygon
#endif

!=======================================================================
!
!    query_triangle ( nside, v1, v2, v3, listpix, nlist, nest, inclusive)
!    --------------
!     nside       = resolution parameter (a power of 2)
!     v1, v2, v3  = 3D vector location of the 3 triangle vertices
!     list_pix    = list of pixel lying in the triangle
!     nlist       = number of pixels in the list
!     nest  (OPT), :0 by default, the output list is in RING scheme
!                  if set to 1, the output list is in NESTED scheme
!     inclusive (OPT) , :0 by default, only the pixels whose center
!                       lie in the triangle are listed on output
!                  if set to 1, all pixels overlapping the triangle are output
!
!     NB : the dimension of the listpix array is fixed in the calling
!     routine and should be large enough for the specific configuration
!
!
! v1.0, EH, Caltech, Nov-Dec-2001
! v1.1,  Aug-Sep-2002 : added nest and inclusive
! v1.2, 2008-03-28, EH/IAP: check that input vectors are 3D
! 2009-06-17: deals with Nside > 8192
  !=======================================================================
#ifdef DOI8B
  subroutine query_triangle_8 ( nside, v1, v2, v3, listpix, nlist, nest, inclusive)
    integer(kind=i4b), parameter  ::   MKD = I8B
#else
  subroutine query_triangle   ( nside, v1, v2, v3, listpix, nlist, nest, inclusive)
    integer(kind=i4b), parameter  ::   MKD = I4B
#endif
    !=======================================================================
    integer(kind=I4B), intent(in) :: nside
    real(kind=DP),     intent(in),  dimension(1:)  :: v1, v2, v3
    integer(kind=MKD), intent(out), dimension(0:)  :: listpix
    integer(kind=MKD), intent(out)                 :: nlist
    integer(kind=I4B), intent(in), optional        :: nest
    integer(kind=I4B), intent(in), optional        :: inclusive

    integer(kind=MKD) :: npix, ilist !, ip1, ip2, ip3
    integer(kind=MKD) :: list_size, nlost
    integer(kind=I4B) :: iz, irmin, irmax
    integer(kind=I4B) :: n12, n123a, n123b, ndom
    integer(kind=I4B), dimension(:),   allocatable  :: listir
    logical(kind=LGT) :: test1, test2, test3
    logical(kind=LGT) :: test1a, test1b, test2a, test2b, test3a, test3b
!!!    real(kind=DP) :: dth1, dth2
    real(kind=DP) :: determ, sdet
    real(kind=DP) :: zmax, zmin, z1max, z1min, z2max, z2min, z3max, z3min
    real(kind=DP) :: z, zz
    real(kind=DP) :: tgth, st
    real(kind=DP) :: offset, sin_off, cos_off
    real(kind=DP), dimension(1:3,1:3) :: vv, vo
    real(kind=DP), dimension(1:3) :: sprod, sto, phi0i, tgthi
    real(kind=DP), dimension(1:3) :: dc
    real(kind=DP), dimension(1:2,1:4) :: dom
    !real(kind=DP), dimension(1:4) :: dom12, dom123a, dom123b
    real(kind=DP), dimension(1:8) :: alldom, adtmp
    real(kind=DP) :: a_i, b_i, phi0, dphiring, dphi, phi0disc, radius_eff, tmp
    integer(kind=I4B) :: idom, nir, ip, nr
    integer(kind=I4B) :: status
    integer(kind=I4B) :: j
    character(len=*), parameter :: code = "QUERY_TRIANGLE"
    logical(LGT)      :: do_inclusive
    integer(kind=I4B)       :: my_nest
    integer(i4b), dimension(1:1) :: longside
    integer(i4b) :: lsp1, lsp2, i
    real(dp), dimension(1:3) :: vcenter, dd
    real(dp), allocatable, dimension(:) :: ztab, dphitab

    !=======================================================================

    list_size = long_size(listpix)

    npix = nside2npix(nside)
    if (npix < 0) then
       print*,"Invalid Nside = ",nside
       call fatal_error
    endif

    do_inclusive = .false.
    if (present(inclusive)) then
       if (inclusive == 1) do_inclusive = .true.
    endif

    my_nest = 0
    if (present(nest)) then
       if (nest == 0 .or. nest == 1) then
          my_nest = nest
       else
          print*,code//"> NEST should be 0 or 1"
          call fatal_error("> program abort ")
       endif
    endif

    if (size(v1) < 3 .or. size(v2) < 3 .or. size(v3) < 3) then
       print*,'Expecting 3d vectors in '//code
       call fatal_error
    endif

    vv(1:3,1) = v1(1:3) / sqrt(dot_product(v1,v1))
    vv(1:3,2) = v2(1:3) / sqrt(dot_product(v2,v2))
    vv(1:3,3) = v3(1:3) / sqrt(dot_product(v3,v3))

    !     --------- allocate memory -------------
    ALLOCATE( listir(0: 4*nside-1), STAT = status)
    if (status /= 0) then
       write(unit=*,fmt="(a)") code//" > can not allocate memory for listir :"
       call fatal_error(" > program abort ")
    endif

!     dth1 = 1.0_dp / (3.0_dp*REAL(nside,kind=dp)**2)
!     dth2 = 2.0_dp / (3.0_dp*REAL(nside,kind=dp))


    ! determ = (vect1 X vect2) . vect3
    ! determines the left(<0)/right(>0) handedness of the triangle
    determ = vv(1,1)*vv(2,2)*vv(3,3) + vv(1,2)*vv(2,3)*vv(3,1) + vv(1,3)*vv(2,1)*vv(3,2) &
         & - vv(3,1)*vv(2,2)*vv(1,3) - vv(3,2)*vv(2,3)*vv(1,1) - vv(3,3)*vv(2,1)*vv(1,2)

    if (abs(determ) < 1.e-20_dp) then
       print*,' ************************************************************'
       print*,' The triangle is degenerated (2 of the vertices are antipodal)'
       print*,' The query can not be performed '
       print*,' ************************************************************'
       call fatal_error
    endif

    !   print*,determ
    sdet = SIGN(1.0_dp,determ) ! = +1 or -1, the sign of determ

    ! scalar product of vertices vectors
    sprod(1) = dot_product(vv(1:3,2),vv(1:3,3))
    sprod(2) = dot_product(vv(1:3,3),vv(1:3,1))
    sprod(3) = dot_product(vv(1:3,1),vv(1:3,2))

    ! vector orthogonal to the great circle containing the vertex doublet
    call vect_prod(vv(1:3,2), vv(1:3,3), vo(1:3,1))
    call vect_prod(vv(1:3,3), vv(1:3,1), vo(1:3,2))
    call vect_prod(vv(1:3,1), vv(1:3,2), vo(1:3,3))

    ! normalize the orthogonal vector
    vo(1:3,1) = vo(1:3,1) /  SQRT(SUM(vo(1:3,1)**2))
    vo(1:3,2) = vo(1:3,2) /  SQRT(SUM(vo(1:3,2)**2))
    vo(1:3,3) = vo(1:3,3) /  SQRT(SUM(vo(1:3,3)**2))
    
    ! test presence of poles in the triangle
    zmax = -1.0_dp
    zmin =  1.0_dp
    test1 = (vo(3,1) * sdet >= 0.0_dp) ! north pole in hemisphere defined by 2-3
    test2 = (vo(3,2) * sdet >= 0.0_dp) ! north pole in hemisphere defined by 1-2
    test3 = (vo(3,3) * sdet >= 0.0_dp) ! north pole in hemisphere defined by 1-3
    if (test1 .and. test2 .and. test3) then
       zmax = 1.0_dp ! north pole in the triangle
    endif
    if ((.not.test1) .and. (.not.test2) .and. (.not.test3)) then
       zmin = -1.0_dp ! south pole in the triangle
    endif

    ! look for northernest and southernest points in the triangle
    !   node(1,2) = vector of norm=1, in the plane defined by (1,2) and with z=0
    test1a = ((vv(3,3) - sprod(1) * vv(3,2)) >= 0.0_dp) ! segment 2-3 : -vector(3) . node(2,3)
    test1b = ((vv(3,2) - sprod(1) * vv(3,3)) >= 0.0_dp) !                vector(2) . node(2,3)
    test2a = ((vv(3,3) - sprod(2) * vv(3,1)) >= 0.0_dp) ! segment 1-3 : -vector(3) . node(1,3)
    test2b = ((vv(3,1) - sprod(2) * vv(3,3)) >= 0.0_dp) !                vector(1) . node(1,3)
    test3a = ((vv(3,2) - sprod(3) * vv(3,1)) >= 0.0_dp) ! segment 1-2 : -vector(2) . node(1,2)
    test3b = ((vv(3,1) - sprod(3) * vv(3,2)) >= 0.0_dp) !                vector(1) . node(1,2)

    ! sin of theta for orthogonal vector
    sto(1:3) = SQRT( (1.0_dp-vo(3,1:3))*(1.0_dp+vo(3,1:3)) )

    ! for each segment (=side of the triangle) the extrema are either
    ! - the 2 vertices
    ! - one of the vertices and a point within the segment

    ! segment 2-3
    z1max = vv(3,2)
    z1min = vv(3,3)
    if ( test1a .EQV. test1b ) then
       zz = sto(1)
       if (z1min+z1max >= 0.0_dp) then
          z1max =  zz
       else
          z1min = -zz
       endif
    endif

    ! segment 1-3
    z2max = vv(3,3)
    z2min = vv(3,1)
    if ( test2a .EQV. test2b ) then
       zz = sto(2)
       if (z2min+z2max >= 0.0_dp) then
          z2max =  zz
       else
          z2min = -zz
       endif
    endif

    ! segment 1-2
    z3max = vv(3,1)
    z3min = vv(3,2)
    if ( test3a .EQV. test3b ) then
       zz = sto(3)
       if (z3min+z3max >= 0.0_dp) then
          z3max =  zz
       else
          z3min = -zz
       endif
    endif

    zmax = MAX(z1max, z2max, z3max, zmax)
    zmin = MIN(z1min, z2min, z3min, zmin)

    offset = 0.0_dp
    sin_off = 0.0_dp
    if (do_inclusive) then
       ! if we are inclusive, move the upper point up, and the lower point down, by a fudge offset
       offset = fudge_query_radius(nside)
       sin_off = sin(offset)
       cos_off = cos(offset)
       zmax = min( 1.0_dp, cos_off * zmax + sin_off * sqrt(1.0_dp - zmax**2))  !cos(theta_zmax-offset)
       zmin = max(-1.0_dp, cos_off * zmin - sin_off * sqrt(1.0_dp - zmin**2))  !cos(theta_zmin+offset)

       ! find one small circle containing all points, 
       ! increased by fudge offset in inclusive case
       longside = minloc(sprod)      ! l   in {1,2,3}: longest leg
       lsp1 = mod(longside(1)  , 3) + 1 ! l+1 in {2,3,1}
       lsp2 = mod(longside(1)+1, 3) + 1 ! l+2 in {3,1,2}
       vcenter = 0.5_dp*(vv(1:3,lsp1) + vv(1:3,lsp2)) ! mid point of longest side
       vcenter = vcenter/sqrt(sum(vcenter**2))
       do i=1,3 
          call angdist(vcenter, vv(1:3,i), dd(i))
       enddo
       radius_eff = maxval(dd) + offset
!        print*,'radius_eff, offset',radius_eff, offset
!        print*,'vcenter',vcenter
    endif


    ! northernmost and sourthernmost ring number
    irmin = ring_num(nside, zmax)
    irmax = ring_num(nside, zmin)

    ilist = -1

    ! -------- loop on the rings -------------------------

    tgthi(1:3) = -1.0e30_dp * vo(3,1:3)
    phi0i(1:3) =  0.0_dp
    do j=1,3
       if (sto(j) > 1.0e-10_dp) then
          tgthi(j) = -vo(3,j) / sto(j)  ! -cotan(theta_orth)
          phi0i(j) = ATAN2(vo(2,j),vo(1,j))
       endif
    enddo
    !   print*,tgthi,phi0i

    ! the triangle boundaries are geodesics : intersection of the sphere with plans going thru (0,0,0)
    ! if we are inclusive, the boundaries are the intersecion of the sphere with plans pushed outward
    ! by sin(offset)

    nr = irmax - irmin + 1
    allocate(ztab(1:nr))
    do iz = irmin, irmax
       ztab(iz-irmin+1) = ring2z(nside, iz)
    enddo
    if (do_inclusive) then
       allocate(dphitab(1:nr))
       call discphirange_at_z(vcenter, radius_eff, ztab, nr, dphitab, phi0disc)
    endif

    do iz = irmin, irmax
       z = ztab(iz-irmin+1) ! z of ring being considered

       ! computes the 3 intervals described by the 3 great circles
       st = SQRT((1.0_dp - z)*(1.0_dp + z))
       tgth = z / st ! cotan(theta_ring)
       dc(1:3)  = tgthi(1:3) * tgth - sdet * sin_off / ((sto(1:3)+1.e-30_dp) * st) ! sto is slightly offset to avoid division by 0

       do j=1,3
          if (dc(j)*sdet <= -1.0_dp) then  ! the whole iso-latitude ring is on the right side of the great circle
             dom(1, j) = 0.0_dp
             dom(2, j) = twopi
          else if (dc(j)*sdet >= 1.0_dp) then ! all on the wrong side
             dom(1, j) = -1.000001_dp * j
             dom(2, j) = -1.0_dp * j
          else ! some is good, some is bad
             tmp = acos(dc(j)) * sdet
             dom(1, j) = modulo( phi0i(j) - tmp , twopi)
             dom(2, j) = modulo( phi0i(j) + tmp , twopi)
          endif
       enddo
       dom(1:2,4) = (/ -2.000002_dp, -2.0_dp /)
       if (do_inclusive) then
          dphi = dphitab(iz-irmin+1)
          if (dphi >= 0.0_dp) then
             dom(1,4) =  modulo (phi0disc - dphi + twopi, twopi)
             dom(2,4) =  modulo (phi0disc + dphi + twopi, twopi)
          endif
       endif

       call process_intervals(dom(1:2,2),dom(1:2,1), alldom, ndom)
       if (ndom == 0) goto 20
       adtmp(1:2*ndom) = alldom(1:2*ndom) ! to avoid using same variable in input and output of process_intervals
       call process_intervals(dom(1:2,3), adtmp(1:2*ndom), alldom, ndom)
       if (ndom == 0) goto 20
       if (do_inclusive) then
          adtmp(1:2*ndom) = alldom(1:2*ndom)
          call process_intervals(dom(1:2,4), adtmp(1:2*ndom), alldom, ndom)
          if (ndom == 0) goto 20
       endif
       do idom=0,ndom-1
          a_i = alldom(2*idom+1)
          b_i = alldom(2*idom+2)
          phi0     = (a_i + b_i) * 0.5_dp
          dphiring = (b_i - a_i) * 0.5_dp
          if (dphiring < 0.0_dp) then
             phi0     = phi0 + pi
             dphiring = dphiring + pi
          endif

          !        ------- finds pixels in the triangle on that ring ---------
          call in_ring(nside, iz, phi0, dphiring, listir, nir, nest=my_nest)

          !        ----------- merge pixel lists -----------
          nlost = ilist + nir + 1_MKD - list_size
          if ( nlost > 0 ) then
             print*,code//"> listpix is too short, it will be truncated at ",nir
             print*,"                         pixels lost : ", nlost
             print*, list_size
             nir = nir - nlost
          endif
          do ip = 0, nir-1
             ilist = ilist + 1
             listpix(ilist) = listir(ip)
          enddo
       enddo
20     continue
    enddo !-----------------------------------------

    !     ------ total number of pixel in the disc --------
    nlist = ilist + 1_MKD

    !     ------- deallocate memory and exit ------
    DEALLOCATE(listir)
    deallocate(ztab)
    if (allocated(dphitab)) deallocate(dphitab)

    return
#ifdef DOI8B
  end subroutine query_triangle_8
#else
  end subroutine query_triangle
#endif

!=======================================================================
!     gives the pixel number ipix (NESTED)
!     corresponding to ix, iy and face_num
!
!     Benjamin D. Wandelt 13/10/97
!     using code from HEALPIX toolkit by K.Gorski and E. Hivon
!     2009-06-15: deals with Nside > 8192
!     2012-03-02: test validity of ix_in and iy_in instead of undefined ix and iy
!=======================================================================
#ifdef DOI8B
  subroutine xy2pix_nest_8(nside, ix_in, iy_in, face_num, ipix)
    integer(kind=i4b), parameter  ::   MKD = I8B
#else
  subroutine xy2pix_nest(nside, ix_in, iy_in, face_num, ipix)
    integer(kind=i4b), parameter  ::   MKD = I4B
#endif
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) ::  nside, ix_in, iy_in, face_num
    INTEGER(KIND=MKD), INTENT(OUT) :: ipix
    INTEGER(KIND=I4B) ::  ix, iy, ix_low, ix_hi, iy_low, iy_hi, i, ismax
    integer(kind=MKD) :: ipf, scale, scale_factor
    character(len=*), parameter :: code = "xy2pix_nest"

    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) call fatal_error(code//"> nside out of range")
    if (ix_in<0 .or. ix_in>(nside-1)) call fatal_error(code//"> ix out of range")
    if (iy_in<0 .or. iy_in>(nside-1)) call fatal_error(code//"> iy out of range")
    if (x2pix1(127) <= 0) call mk_xy2pix1()

    ix = ix_in
    iy = iy_in
    if (nside <= ns_max4) then 
       ix_low = iand(ix, 127)
       iy_low = iand(iy, 127)
       ipf =     x2pix1(ix_low) + y2pix1(iy_low) &
            & + (x2pix1(ix/128) + y2pix1(iy/128)) * 16384
    else
       scale = 1_MKD
       scale_factor = 16384_MKD ! 128*128
       ipf = 0_MKD
       ismax = 1 ! for nside in [2^14, 2^20]
       if (nside >  1048576 ) ismax = 3
       do i=0, ismax
          ix_low = iand(ix, 127) ! last 7 bits
          iy_low = iand(iy, 127) ! last 7 bits
          ipf = ipf + (x2pix1(ix_low)+y2pix1(iy_low)) * scale
          scale = scale * scale_factor
          ix  =     ix / 128 ! truncate out last 7 bits
          iy  =     iy / 128
       enddo
       ipf =  ipf + (x2pix1(ix)+y2pix1(iy)) * scale
    endif
    ipix = ipf + face_num* int(nside,MKD) * nside    ! in {0, 12*nside**2 - 1}
    return

#ifdef DOI8B
  end subroutine xy2pix_nest_8
#else
  end subroutine xy2pix_nest
#endif
!=======================================================================
!  pix2xy_nest
!     gives the x, y coords in a face from pixel number within the face (NESTED)
!
!     Benjamin D. Wandelt 13/10/97
!
!     using code from HEALPIX toolkit by K.Gorski and E. Hivon
!     2009-06-15: deals with Nside > 8192
!     2012-03-02: test validity of ipf_in instead of undefined ipf
!                 define ipf as MKD
!     2012-08-27:  corrected bug on (ix,iy) for Nside > 8192 (MARK)
!=======================================================================
#ifdef DOI8B
  subroutine pix2xy_nest_8(nside, ipf_in, ix, iy)
    integer(kind=i4b), parameter  ::   MKD = I8B
#else
  subroutine pix2xy_nest  (nside, ipf_in, ix, iy)
    integer(kind=i4b), parameter  ::   MKD = I4B
#endif
    INTEGER(KIND=I4B), INTENT(IN)  :: nside
    INTEGER(KIND=MKD), INTENT(IN)  :: ipf_in
    INTEGER(KIND=I4B), INTENT(OUT) :: ix, iy

    integer(kind=MKD) :: ipf
    INTEGER(KIND=I4B) ::  ip_low, ip_trunc, ip_med, ip_hi, scale, i, ismax
    character(len=*), parameter :: code = "pix2xy_nest"

    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) call fatal_error(code//"> nside out of range")
    if (ipf_in<0 .or. ipf_in>nside*nside-1) &
         &     call fatal_error(code//"> ipix out of range")
    if (pix2x(1023) <= 0) call mk_pix2xy()

    ipf = ipf_in
    if (nside <= ns_max4) then
       ip_low = iand(ipf,1023_MKD)   ! content of the last 10 bits
       ip_trunc =    ipf/1024        ! truncation of the last 10 bits
       ip_med = iand(ip_trunc,1023)  ! content of the next 10 bits
       ip_hi  =      ip_trunc/1024   ! content of the high weight 10 bits

       ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
       iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)
    else
       ix = 0
       iy = 0
       scale = 1
       ismax = 4
       do i=0, ismax
          ip_low = iand(ipf,1023_MKD)
          ix = ix + scale * pix2x(ip_low)
          iy = iy + scale * pix2y(ip_low) ! corrected 2012-08-27
          scale = scale * 32
          ipf   = ipf/1024
       enddo
       ix = ix + scale * pix2x(ipf)
       iy = iy + scale * pix2y(ipf) ! corrected 2012-08-27
    endif

    return

#ifdef DOI8B
  end subroutine pix2xy_nest_8
#else
  end subroutine pix2xy_nest
#endif



!====================================================================
! The following is a routine which finds the 7 or 8 neighbours of
! any pixel in the nested scheme of the HEALPIX pixelisation.
!====================================================================
!  neighbours_nest
!
!   Returns list n(8) of neighbours of pixel ipix (in NESTED scheme)
!   the neighbours are ordered in the following way:
!   First pixel is the one to the south (the one west of the south
! direction is taken
! for the pixels which don't have a southern neighbour). From
! then on the neighbours are ordered in the clockwise direction
! about the pixel with number ipix.
!
!   nneigh is the number of neighbours (mostly 8, 8 pixels have 7 neighbours)
!
!   Benjamin D. Wandelt October 1997
!   Added to pix_tools in March 1999
!   added 'return' for case nside=1, EH, Oct 2005
!   corrected bugs in case nside=1 and ipix=7, 9 or 11, EH, June 2006
!   2009-06-16: deals with Nside > 8192
!====================================================================
#ifdef DOI8B
  subroutine neighbours_nest_8(nside, ipix, n, nneigh)
    use bit_manipulation
    integer(kind=i4b), parameter  ::   MKD = I8B
#else
  subroutine neighbours_nest(nside, ipix, n, nneigh)
    use bit_manipulation
    integer(kind=i4b), parameter  ::   MKD = I4B
#endif
    !====================================================================
    integer(kind=i4b), intent(in)::  nside
    integer(kind=MKD), intent(in)::  ipix
    integer(kind=MKD), intent(out), dimension(1:):: n
    integer(kind=i4b), intent(out):: nneigh

    integer(kind=i4b) :: ix,ixm,ixp,iy,iym,iyp,ixo,iyo
    integer(kind=i4b) :: face_num,other_face, ipix4
    integer(kind=i4b) :: ia,ib,ibp,ibm,ib2,icase
    integer(kind=MKD) :: npix,ipf,ipo
    integer(kind=MKD) :: local_magic1,local_magic2,nsidesq
    integer(kind=i4b), dimension(1:8) :: n4
    character(len=*), parameter :: code = "neighbours_nest"

!     integer(kind=i4b), intrinsic :: IAND

    !--------------------------------------------------------------------
#ifdef DOI8B
    if (nside <= ns_max4) then ! use faster 32-bit routine whenever possible
       ipix4 = int(ipix, kind=i4b)
       call neighbours_nest(nside, ipix4, n4, nneigh)
       n(1:nneigh) = int(n4(1:nneigh), kind=MKD)
       return
    endif
    if (nside <1 .or. nside > ns_max ) call fatal_error(code//"> nside out of range")
#else
    if (nside <1 .or. nside > ns_max4) call fatal_error(code//"> nside out of range")
#endif
    npix = nside2npix(nside) ! total number of points
    nsidesq = npix / 12
    if (ipix <0 .or. ipix>npix-1) call fatal_error(code//"> ipix out of range")

    ! quick and dirty hack for Nside=1

    if (nside == 1) then
       nneigh = 6
       if (ipix==0 ) n(1:6) = (/ 8, 4, 3, 2, 1, 5 /)
       if (ipix==1 ) n(1:6) = (/ 9, 5, 0, 3, 2, 6 /)
       if (ipix==2 ) n(1:6) = (/10, 6, 1, 0, 3, 7 /)
       if (ipix==3 ) n(1:6) = (/11, 7, 2, 1, 0, 4 /)
       if (ipix==4 ) n(1:6) = (/11, 7, 3, 0, 5, 8 /)
       if (ipix==5 ) n(1:6) = (/ 8, 4, 0, 1, 6, 9 /)
       if (ipix==6 ) n(1:6) = (/ 9, 5, 1, 2, 7,10 /)
       if (ipix==7 ) n(1:6) = (/10, 6, 2, 3, 4,11 /)
       if (ipix==8 ) n(1:6) = (/10,11, 4, 0, 5, 9 /)
       if (ipix==9 ) n(1:6) = (/11, 8, 5, 1, 6,10 /)
       if (ipix==10) n(1:6) = (/ 8, 9, 6, 2, 7,11 /)
       if (ipix==11) n(1:6) = (/ 9,10, 7, 3, 4, 8 /)
       return
    endif

    !     initiates array for (x,y)-> pixel number -> (x,y) mapping
    if (x2pix1(127) <= 0) call mk_xy2pix1()

    local_magic1=(nsidesq-1)/3
    local_magic2=2*local_magic1
    face_num=ipix/nsidesq

    ipf=modulo(ipix,nsidesq)   !Pixel number in face

    call pix2xy_nest(nside,ipf,ix,iy)
    ixm=ix-1
    ixp=ix+1
    iym=iy-1
    iyp=iy+1

    nneigh=8                  !Except in special cases below

    !     Exclude corners
    if(ipf==local_magic2)     then !WestCorner
       icase=5
       goto 100
    endif
    if(ipf==(nsidesq-1)) then !NorthCorner
       icase=6
       goto 100
    endif
    if(ipf==0)           then !SouthCorner
       icase=7
       goto 100
    endif
    if(ipf==local_magic1)     then !EastCorner
       icase=8
       goto 100
    endif

    !     Detect edges
    if(IAND(ipf,local_magic1)==local_magic1) then !NorthEast
       icase=1
       goto 100
    endif
    if(IAND(ipf,local_magic1)==0)      then !SouthWest
       icase=2
       goto 100
    endif
    if(IAND(ipf,local_magic2)==local_magic2) then !NorthWest
       icase=3
       goto 100
    endif
    if(IAND(ipf,local_magic2)==0)      then !SouthEast
       icase=4
       goto 100
    endif

    !     Inside a face
    call xy2pix_nest(nside, ixm, iym, face_num, n(1))
    call xy2pix_nest(nside, ixm, iy , face_num, n(2))
    call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
    call xy2pix_nest(nside, ix , iyp, face_num, n(4))
    call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
    call xy2pix_nest(nside, ixp, iy , face_num, n(6))
    call xy2pix_nest(nside, ixp, iym, face_num, n(7))
    call xy2pix_nest(nside, ix , iym, face_num, n(8))
    return

100 continue

    ia= face_num/4            !in {0,2}
    ib= modulo(face_num,4)       !in {0,3}
    ibp=modulo(ib+1,4)
    ibm=modulo(ib+4-1,4)
    ib2=modulo(ib+2,4)

    if(ia==0) then          !North Pole region
       select case(icase)
       case(1)              !NorthEast edge
          other_face=0+ibp
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ixm, iym, face_num, n(1))
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          ipo=modulo(swapLSBMSB(ipf),nsidesq)    !East-West flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo+1 , iyo, other_face, n(5))
          n(6)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo-1, iyo, other_face, n(7))
       case(2)              !SouthWest edge
          other_face=4+ib
          ipo=modulo(invLSB(ipf),nsidesq)        !SW-NE flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo, iyo-1, other_face, n(1))
          n(2)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo, iyo+1, other_face, n(3))
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          call xy2pix_nest(nside, ixp, iym, face_num, n(7))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
          call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
       case(3)              !NorthWest edge
          other_face=0+ibm
          ipo=modulo(swapLSBMSB(ipf),nsidesq)    !East-West flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo, iyo-1, other_face, n(3))
          n(4)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo, iyo+1, other_face, n(5))
          call xy2pix_nest(nside, ixm, iym, face_num, n(1))
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ixp, iym, face_num, n(7))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
       case(4)              !SouthEast edge
          other_face=4+ibp
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
          ipo=modulo(invMSB(ipf),nsidesq) !SE-NW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo+1, iyo, other_face, n(7))
          n(8)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo-1, iyo, other_face, n(1))
       case(5)              !West corner
          nneigh=7
          other_face=4+ib
          n(2)=other_face*nsidesq+nsidesq-1
          n(1)=n(2)-2
          other_face=0+ibm
          n(3)=other_face*nsidesq+local_magic1
          n(4)=n(3)+2
          n(5)=ipix+1
          n(6)=ipix-1
          n(7)=ipix-2
       case(6)              !North corner
          n(1)=ipix-3
          n(2)=ipix-1
          n(8)=ipix-2
          other_face=0+ibm
          n(4)=other_face*nsidesq+nsidesq-1
          n(3)=n(4)-2
          other_face=0+ib2
          n(5)=other_face*nsidesq+nsidesq-1
          other_face=0+ibp
          n(6)=other_face*nsidesq+nsidesq-1
          n(7)=n(6)-1
       case(7)              !South corner
          other_face=8+ib
          n(1)=other_face*nsidesq+nsidesq-1
          other_face=4+ib
          n(2)=other_face*nsidesq+local_magic1
          n(3)=n(2)+2
          n(4)=ipix+2
          n(5)=ipix+3
          n(6)=ipix+1
          other_face=4+ibp
          n(8)=other_face*nsidesq+local_magic2
          n(7)=n(8)+1
       case(8)              !East corner
          nneigh=7
          n(2)=ipix-1
          n(3)=ipix+1
          n(4)=ipix+2
          other_face=0+ibp
          n(6)=other_face*nsidesq+local_magic2
          n(5)=n(6)+1
          other_face=4+ibp
          n(7)=other_face*nsidesq+nsidesq-1
          n(1)=n(7)-1
       end select ! north

    elseif(ia==1) then      !Equatorial region
       select case(icase)
       case(1)              !NorthEast edge
          other_face=0+ib
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ixm, iym, face_num, n(1))
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          ipo=modulo(invLSB(ipf),nsidesq)    !NE-SW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo , iyo+1, other_face, n(5))
          n(6)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo, iyo-1, other_face, n(7))
       case(2)              !SouthWest edge
          other_face=8+ibm
          ipo=modulo(invLSB(ipf),nsidesq)        !SW-NE flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo, iyo-1, other_face, n(1))
          n(2)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo, iyo+1, other_face, n(3))
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          call xy2pix_nest(nside, ixp, iym, face_num, n(7))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
          call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
       case(3)              !NorthWest edge
          other_face=0+ibm
          ipo=modulo(invMSB(ipf),nsidesq)    !NW-SE flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo-1, iyo, other_face, n(3))
          n(4)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo+1, iyo, other_face, n(5))
          call xy2pix_nest(nside, ixm, iym, face_num, n(1))
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ixp, iym, face_num, n(7))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
       case(4)              !SouthEast edge
          other_face=8+ib
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
          ipo=modulo(invMSB(ipf),nsidesq) !SE-NW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo+1, iyo, other_face, n(7))
          n(8)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo-1, iyo, other_face, n(1))
       case(5)              !West corner
          other_face=8+ibm
          n(2)=other_face*nsidesq+nsidesq-1
          n(1)=n(2)-2
          other_face=4+ibm
          n(3)=other_face*nsidesq+local_magic1
          other_face=0+ibm
          n(4)=other_face*nsidesq
          n(5)=n(4)+1
          n(6)=ipix+1
          n(7)=ipix-1
          n(8)=ipix-2
       case(6)              !North corner
          nneigh=7
          n(1)=ipix-3
          n(2)=ipix-1
          other_face=0+ibm
          n(4)=other_face*nsidesq+local_magic1
          n(3)=n(4)-1
          other_face=0+ib
          n(5)=other_face*nsidesq+local_magic2
          n(6)=n(5)-2
          n(7)=ipix-2
       case(7)              !South corner
          nneigh=7
          other_face=8+ibm
          n(1)=other_face*nsidesq+local_magic1
          n(2)=n(1)+2
          n(3)=ipix+2
          n(4)=ipix+3
          n(5)=ipix+1
          other_face=8+ib
          n(7)=other_face*nsidesq+local_magic2
          n(6)=n(7)+1
       case(8)              !East corner
          other_face=8+ib
          n(8)=other_face*nsidesq+nsidesq-1
          n(1)=n(8)-1
          n(2)=ipix-1
          n(3)=ipix+1
          n(4)=ipix+2
          other_face=0+ib
          n(6)=other_face*nsidesq
          n(5)=n(6)+2
          other_face=4+ibp
          n(7)=other_face*nsidesq+local_magic2
       end select ! equator
    else                    !South Pole region
       select case(icase)
       case(1)              !NorthEast edge
          other_face=4+ibp
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ixm, iym, face_num, n(1))
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          ipo=modulo(invLSB(ipf),nsidesq)    !NE-SW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo , iyo+1, other_face, n(5))
          n(6)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo, iyo-1, other_face, n(7))
       case(2)              !SouthWest edge
          other_face=8+ibm
          ipo=modulo(swapLSBMSB(ipf),nsidesq)        !W-E flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo-1, iyo, other_face, n(1))
          n(2)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo+1, iyo, other_face, n(3))
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          call xy2pix_nest(nside, ixp, iym, face_num, n(7))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
          call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
       case(3)              !NorthWest edge
          other_face=4+ib
          ipo=modulo(invMSB(ipf),nsidesq)    !NW-SE flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo-1, iyo, other_face, n(3))
          n(4)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo+1, iyo, other_face, n(5))
          call xy2pix_nest(nside, ixm, iym, face_num, n(1))
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ixp, iym, face_num, n(7))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
       case(4)              !SouthEast edge
          other_face=8+ibp
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
          ipo=modulo(swapLSBMSB(ipf),nsidesq) !E-W flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo, iyo+1, other_face, n(7))
          n(8)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo, iyo-1, other_face, n(1))
       case(5)              !West corner
          nneigh=7
          other_face=8+ibm
          n(2)=other_face*nsidesq+local_magic1
          n(1)=n(2)-1
          other_face=4+ib
          n(3)=other_face*nsidesq
          n(4)=n(3)+1
          n(5)=ipix+1
          n(6)=ipix-1
          n(7)=ipix-2
       case(6)              !North corner
          n(1)=ipix-3
          n(2)=ipix-1
          other_face=4+ib
          n(4)=other_face*nsidesq+local_magic1
          n(3)=n(4)-1
          other_face=0+ib
          n(5)=other_face*nsidesq
          other_face=4+ibp
          n(6)=other_face*nsidesq+local_magic2
          n(7)=n(6)-2
          n(8)=ipix-2
       case(7)              !South corner
          other_face=8+ib2
          n(1)=other_face*nsidesq
          other_face=8+ibm
          n(2)=other_face*nsidesq
          n(3)=n(2)+1
          n(4)=ipix+2
          n(5)=ipix+3
          n(6)=ipix+1
          other_face=8+ibp
          n(8)=other_face*nsidesq
          n(7)=n(8)+2
       case(8)              !East corner
          nneigh=7
          other_face=8+ibp
          n(7)=other_face*nsidesq+local_magic2
          n(1)=n(7)-2
          n(2)=ipix-1
          n(3)=ipix+1
          n(4)=ipix+2
          other_face=4+ibp
          n(6)=other_face*nsidesq
          n(5)=n(6)+2
       end select ! south
    endif

    return

#ifdef DOI8B
  end subroutine neighbours_nest_8
#else
  end subroutine neighbours_nest
#endif
!====================================================================
!  next_in_line_nest
!
!   given nside and a NESTED pixel number ipix, returns in inext
!  the pixel that lies on the East side (and the same latitude) as ipix
!
!   Hacked by EH from BDW's neighbours_nest, 2001-12-18
!   Hacked for Nside=1 by EH, 2004-05-28
  !====================================================================
#ifdef DOI8B
  subroutine next_in_line_nest_8(nside, ipix, inext)
    use bit_manipulation
    integer(i4b), parameter :: MKD = I8B
#else
  subroutine next_in_line_nest  (nside, ipix, inext)
    use bit_manipulation
    integer(i4b), parameter :: MKD = I4B
#endif
    !====================================================================
    integer(kind=i4b), intent(in)::nside
    integer(kind=MKD), intent(in)::  ipix
    integer(kind=MKD), intent(out):: inext

    integer(kind=i4b) :: ix,ixp,iy,iym,ixo,iyo
    integer(kind=i4b) :: face_num,other_face
    integer(kind=i4b) :: ia,ib,ibp,ibm,ib2,icase
    integer(kind=i4b) :: ipix4, inext4
    integer(kind=MKD) :: npix,ipf,ipo
    integer(kind=MKD) :: local_magic1,local_magic2,nsidesq
    character(len=*), parameter :: code = "next_in_line_nest"

    !--------------------------------------------------------------------
#ifdef DOI8B
    if (nside <= ns_max4) then ! use faster 32-bit routine whenever possible
       ipix4 = int(ipix, kind=i4b)
       call next_in_line_nest(nside, ipix4, inext4)
       inext = int(inext4, kind=MKD)
       return
    endif
    if (nside <1 .or. nside > ns_max ) call fatal_error(code//"> nside out of range")
#else
    if (nside <1 .or. nside > ns_max4) call fatal_error(code//"> nside out of range")
#endif
    npix = nside2npix(nside) ! total number of points
    nsidesq = npix / 12
    if (ipix <0 .or. ipix>npix-1) call fatal_error(code//"> ipix out of range")

    ! quick and dirty hack for Nside=1
    if (nside == 1) then
       inext = ipix + 1
       if (ipix == 3)  inext = 0_MKD
       if (ipix == 7)  inext = 4_MKD
       if (ipix == 11) inext = 8_MKD
       return
    endif
    !     initiates array for (x,y)-> pixel number -> (x,y) mapping
    if (x2pix1(127) <= 0) call mk_xy2pix1()

    local_magic1=(nsidesq-1)/3
    local_magic2=2*local_magic1
    face_num=ipix/nsidesq

    ipf=modulo(ipix,nsidesq)   !Pixel number in face

    call pix2xy_nest(nside,ipf,ix,iy)
    ixp=ix+1
    iym=iy-1

    !     Exclude corners
    if(ipf==local_magic2)     then !WestCorner
       inext = ipix - 1
       return
    endif
    if(ipf==(nsidesq-1)) then !NorthCorner
       icase=6
       goto 100
    endif
    if(ipf==0)           then !SouthCorner
       icase=7
       goto 100
    endif
    if(ipf==local_magic1)     then !EastCorner
       icase=8
       goto 100
    endif

    !     Detect edges
    if(IAND(ipf,local_magic1)==local_magic1) then !NorthEast
       icase=1
       goto 100
    endif
    if(IAND(ipf,local_magic2)==0)      then !SouthEast
       icase=4
       goto 100
    endif

    !     Inside a face
    call xy2pix_nest(nside, ixp, iym, face_num, inext)
    return

100 continue

    ia= face_num/4            !in {0,2}
    ib= modulo(face_num,4)       !in {0,3}
    ibp=modulo(ib+1,4)
    ibm=modulo(ib+4-1,4)
    ib2=modulo(ib+2,4)

    if(ia==0) then          !North Pole region
       select case(icase)
       case(1)              !NorthEast edge
          other_face=0+ibp
          ipo=modulo(swapLSBMSB(ipf),nsidesq)    !East-West flip
          inext = other_face*nsidesq+ipo         ! (6)
       case(4)              !SouthEast edge
          other_face=4+ibp
          ipo=modulo(invMSB(ipf),nsidesq) !SE-NW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo+1, iyo, other_face, inext)
       case(6)              !North corner
          other_face=0+ibp
          inext=other_face*nsidesq+nsidesq-1
       case(7)              !South corner
          other_face=4+ibp
          inext=other_face*nsidesq+local_magic2+1
       case(8)              !East corner
          other_face=0+ibp
          inext=other_face*nsidesq+local_magic2
       end select ! north

    elseif(ia==1) then      !Equatorial region
       select case(icase)
       case(1)              !NorthEast edge
          other_face=0+ib
          ipo=modulo(invLSB(ipf),nsidesq)    !NE-SW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo, iyo-1, other_face, inext)
       case(4)              !SouthEast edge
          other_face=8+ib
          ipo=modulo(invMSB(ipf),nsidesq) !SE-NW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo+1, iyo, other_face, inext)
       case(6)              !North corner
          other_face=0+ib
          inext=other_face*nsidesq+local_magic2-2
       case(7)              !South corner
          other_face=8+ib
          inext=other_face*nsidesq+local_magic2+1
       case(8)              !East corner
          other_face=4+ibp
          inext=other_face*nsidesq+local_magic2
       end select ! equator
    else                    !South Pole region
       select case(icase)
       case(1)              !NorthEast edge
          other_face=4+ibp
          ipo=modulo(invLSB(ipf),nsidesq)    !NE-SW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo, iyo-1, other_face, inext)
       case(4)              !SouthEast edge
          other_face=8+ibp
          ipo=modulo(swapLSBMSB(ipf),nsidesq) !E-W flip
          inext = other_face*nsidesq+ipo   ! (8)
       case(6)              !North corner
          other_face=4+ibp
          inext=other_face*nsidesq+local_magic2 -2
       case(7)              !South corner
          other_face=8+ibp
          inext=other_face*nsidesq
       case(8)              !East corner
          other_face=8+ibp
          inext=other_face*nsidesq+local_magic2
       end select ! south
    endif

    return
#ifdef DOI8B
  end subroutine next_in_line_nest_8
#else
  end subroutine next_in_line_nest
#endif


    !=======================================================================
    ! template_pixel_ring
    !    Returns template-pixel index corresponding to each (RING ordered) pixel provided
    !    The reflexion (E-W, or N-W or both) to apply to the pixel to match the
    !    template is optionally returned
    !-------------------------------------------------------------
    !   For a given resolution Nside, there are 12*Nside*Nside pixels on the
    !   sphere, whose shape belong to one of the (Nside+6)*Nside/4 templates.
    !   (2 in the case Nside=1).
    !   Any pixel can be matched to a single of these templates by a combination of
    !    - a rotation around the polar axis with 
    !    - reflexion(s) around a meridian and/or the equator.
    !   The reflexion number returned by this routine for each pixel is
    !   Reflexion         Operation to do to match template
    !     0                 rotation only
    !     1                 rotation + East-West swap
    !     2                 rotation + North-South swap
    !     3                 rotation + East-West + North-South swaps
    !   
    !   The template pixels are all located in the Northern Hemisphere, or on the
    !   Equator.
    !   They are chosen to have
    !       z >= 2/3,      0< phi <= Pi/2               [Nside*(Nside+2)/4]
    !       2/3 > z >= 0,  phi=0, or phi=Pi/(4*Nside)   [Nside]
    !   They are numbered from 0, starting at the North Pole, with the index
    !   increasing in phi, and then increasing for decreasing z.
  !=======================================================================
#ifdef DOI8B
  subroutine template_pixel_ring_8(nside, pixel, templ8, reflexion)
    integer(i4b), parameter :: MKD = I8B
#else
  subroutine template_pixel_ring  (nside, pixel, templ8, reflexion)
    integer(i4b), parameter :: MKD = I4B
#endif
    !-------------------------------------------------------------
    character(len=*),           parameter :: code = "Template_pixel_ring"
    integer(I4B),    intent(IN)           :: nside
    integer(MKD),    intent(IN)           :: pixel
    integer(MKD),    intent(OUT)          :: templ8
    integer(I4B),    intent(OUT), optional:: reflexion
    ! local variables
    integer(MKD) :: npix, ncap, n0, ip
    integer(MKD) :: iring ! I4B would be enough, but more convenient for iring*iring
    integer(I4B) :: ns4, ifi, ifd, iff, refl
    real(DP)     :: fip
    !-------------------------------------------------------------

    npix = nside2npix(nside)
    if (npix < 0) then
       print*,'Invalid Nside = ',nside
       call fatal_error(code//' Abort')
    endif

    if (pixel < 0 .or. pixel >= npix) then
       print*,'Invalid Pixel = ',pixel
       call fatal_error(code//' Abort')
    endif

    ncap = 2*nside*(nside+1)
    ns4  = 4*nside
!     n0   = max( (nside*(nside+2))/4, 1)
    n0   = (1_MKD + nside * (nside + 2_MKD)) / 4_MKD


    refl = 0 ! no swap

    ! North polar cap
    if (pixel < ncap) then
       ip    = pixel + 1
       fip   = real(int(ip/2, kind=MKD) , kind=dp)
       iring = int(sqrt(ip/2.0_dp - sqrt(fip) )) + 1 ! counted from NORTH pole, starting at 1
       ifi   = MOD( (ip - 1 - 2*iring*(iring-1)), iring) ! in [0,iring-1], upwards
       ifd   = iring - 1 - ifi                           ! in [0,iring-1], downwards
       iff   = min(ifi, ifd)                             ! in [0,(iring-1)/2]
       templ8 = (iring*iring)/4 + iff
       if (ifd < ifi) refl = 1 ! East-West swap


    ! North-Equatorial region
    elseif (pixel < (npix+ns4)/2) then
       templ8 = n0 + (pixel - ncap)/ns4

    ! South-Equatorial region
    elseif (pixel < (npix-ncap)) then
       templ8 = n0 + (npix - 1 - pixel - ncap)/ns4
       refl = 2  ! North-South swap

    ! South polar cap
    else
       ip =  npix - pixel
       fip   = real(int(ip/2, kind=MKD) , kind=dp)
       iring = int(sqrt(ip/2.0_dp - sqrt(fip) )) + 1 ! counted from SOUTH pole, starting at 1
       ifi   = MOD( (2*iring*(iring+1) - ip), iring)    ! in [0,iring-1], upwards
       ifd  = iring - 1 - ifi                           ! in [0,iring-1], downwards
       iff  = min(ifi, ifd)                             ! in [0,(iring-1)/2]
       templ8 = (iring*iring)/4 + iff
       if (ifd < ifi) refl = 1 ! East-West swap
       refl = refl + 2    ! North-South swap
    endif

    if (present(reflexion)) reflexion = refl

    return
#ifdef DOI8B
  end subroutine template_pixel_ring_8
#else
  end subroutine template_pixel_ring
#endif


  !=======================================================================
  ! same_shape_pixels_ring
  !   Returns the list of RING ordered pixels having the same shape as the
  !   template provided
  !=======================================================================
#ifdef DOI8B
  subroutine same_shape_pixels_ring_8(nside, templ8, list, reflexion, nrep)
    integer(i4b), parameter :: MKD = I8B
#else
  subroutine same_shape_pixels_ring  (nside, templ8, list, reflexion, nrep)
    integer(i4b), parameter :: MKD = I4B
#endif
    !-------------------------------------------------------------

    character(len=*), parameter :: code = "same_shape_pixels_ring"
    integer(I4B),    intent(IN)                        :: nside
    integer(MKD),    intent(IN)                        :: templ8
    integer(MKD),    pointer,   dimension(:), optional :: list
    integer(I4B),    pointer,   dimension(:), optional :: reflexion
    integer(MKD),    intent(OUT)            , optional :: nrep

    integer(MKD) :: npix, ncap, ntplt, n0, nrep_in, iring, in0, is0
    integer(i4b) :: ns4 ! = 4*Nside <= 2^30
    integer(i4b) :: i, ifi
    integer(i4b), dimension(0:3) :: ramp4 = (/ 0, 1, 2, 3 /)
    integer(i4b), dimension(0:7) :: alt8  = (/ 0, 1, 0, 1, 0, 1, 0, 1 /)
    integer(i4b), dimension(0:7) :: rr8
    logical(lgt) :: do_list, do_rfx, do_nrep
    integer(i4b) :: ll=1, lr=1
    !-------------------------------------------------------------

    do_list = (present(list))
    do_rfx  = (present(reflexion))
    do_nrep = (present(nrep))

    if (do_rfx .and. .not. do_list) then
       print*,'Error in '//code//'. Can not have Reflexion without pixel list'
       call fatal_error()
    endif

    npix = nside2npix(nside)
    if (npix < 0) then
       print*,'Invalid Nside = ',nside
       call fatal_error(code//' Abort')
    endif

    ntplt= nside2ntemplates(nside)

    if (templ8 < 0 .or. templ8 >= ntplt) then
       print*,'Error on template argument'
       print*,'Nside = ',nside,', Template = ',templ8
       print*,'Template should be in [0, (1+Nside*(Nside+6))/4-1=',ntplt-1,']'
       call fatal_error(code//' Abort')
    endif

    ncap = 2*nside*(nside+1)
    ns4  = 4*nside
!    n0   = max( (nside*(nside+2))/4, 1)
    n0 = (1_MKD + nside * (nside + 2_MKD))/ 4

    ! find number of repetition of template around the sphere
    if (templ8 >= n0) then
       nrep_in = ns4 * 2_MKD ! Non- Equator: nrep = 8*Nside
       if (templ8 == ntplt - 1) nrep_in = ns4 ! Equator : nrep = 4*Nside
    else
       nrep_in  = 16
       iring = int(sqrt(4*templ8+1.0_dp), kind=MKD) ! in [1,Nside]
       ifi   = templ8 - (iring*iring)/4
       if (2*ifi+1 == iring) nrep_in = 8
    endif
    if (do_nrep) nrep = nrep_in
          

    ! create or expand (or leave alone) list array
    if (do_list) then
       if (associated(list)) then
          ll = lbound(list,dim=1)
          if (long_size(list) < nrep_in) then
             deallocate(list)
             allocate(list(ll:ll+nrep_in-1))
          endif
       else
          allocate(list(ll:ll+nrep_in-1))
       endif
       list = -1 ! default value for unused array elements
    endif

    ! create or expand (or leave alone) reflexion array
    if (do_rfx) then
       if (associated(reflexion)) then
          lr = lbound(reflexion,dim=1)
          if (long_size(reflexion) < nrep_in) then
             deallocate(reflexion)
             allocate(reflexion(lr:lr+nrep_in-1))
          endif
       else
          allocate(reflexion(lr:lr+nrep_in-1))
       endif
       reflexion = -1 ! default value for unused array elements
    endif


    ! fill up list and reflexion arrays
    if (templ8 >= n0) then 
       ! Non Polar pixels (most frequent)
    
       in0 = ncap + (templ8 - n0)*ns4

       if (templ8 == ntplt-1) then 
          ! Equator : nrep = 4*Nside
          if (do_list) list      = in0 + (/ (i, i=0, ns4-1) /)
          if (do_rfx)  reflexion = 0
       else 
          ! Non- Equator: nrep = 8*Nside
          is0 = npix - in0 - ns4
          if (do_list) then
             list(ll    :ll+  ns4-1) = in0 + (/ (i, i=0, ns4-1) /)
             list(ll+ns4:ll+2*ns4-1) = is0 + (/ (i, i=0, ns4-1) /)
          endif
          if (do_rfx) then
             reflexion(lr    :lr+  ns4-1) = 0
             reflexion(lr+ns4:lr+2*ns4-1) = 2
          endif
       endif

    else 
       ! Polar pixels
       in0   = 2*iring*(iring-1)
       is0   = npix - in0 - 4*iring

       if (2*ifi+1 == iring) then 
          ! pixel on phi = Pi/4, Nrep = 8
          if (do_list) then
             list(ll  :ll+3) = in0 + ifi + ramp4*iring
             list(ll+4:ll+7) = is0 + ifi + ramp4*iring
          endif
          if (do_rfx) then
             reflexion = (/ 0, 0, 0, 0, 2, 2, 2, 2 /)
          endif
       else 
          ! Nrep = 16
          rr8(0:6:2) = (/ 0, 1, 2, 3 /) * iring + ifi
          rr8(1:7:2) = (/ 1, 2, 3, 4 /) * iring - ifi - 1
          if (do_list) then
             list(ll  :ll+7 ) = in0 + rr8
             list(ll+8:ll+15) = is0 + rr8
          endif
          if (do_rfx) then
             reflexion(lr  :lr+7 ) = alt8
             reflexion(lr+8:lr+15) = alt8 + 2
          endif
       endif
    endif

    return

#ifdef DOI8B
  end subroutine same_shape_pixels_ring_8
#else
  end subroutine same_shape_pixels_ring
#endif


  !=======================================================================
  ! template_pixel_nest
  !    Returns template-pixel index corresponding to each (NESTED ordered) pixel provided
  !    The reflexion (E-W, or N-W or both) to apply to the pixel to match the
  !    template is optionally returned
  !=======================================================================
#ifdef DOI8B
  subroutine template_pixel_nest_8(nside, pixel, templ8, reflexion)
    integer(i4b), parameter :: MKD = I8B
#else
  subroutine template_pixel_nest  (nside, pixel, templ8, reflexion)
    integer(i4b), parameter :: MKD = I4B
#endif
    character(len=*), parameter :: code = "Template_pixel_nest"
    integer(I4B),    intent(IN)           :: nside
    integer(MKD),    intent(IN)           :: pixel
    integer(MKD),    intent(OUT)          :: templ8
    integer(I4B),    intent(OUT), optional:: reflexion
    !
    integer(MKD) :: pix_ring, npix
    
    npix = nside2npix(nside)
    if (npix < 0) then
       print*,'Invalid Nside = ',nside
       call fatal_error(code//' Abort')
    endif

    if (pixel < 0 .or. pixel >= npix) then
       print*,'Invalid Pixel = ',pixel
       call fatal_error(code//' Abort')
    endif

    call nest2ring(nside, pixel, pix_ring)
    call template_pixel_ring(nside, pix_ring, templ8, reflexion)

    return

#ifdef DOI8B
  end subroutine template_pixel_nest_8
#else
  end subroutine template_pixel_nest
#endif

  !=======================================================================
  ! same_shape_pixels_nest  
  !   Returns the list of NESTED ordered pixels having the same shape as the
  !   template provided
  !=======================================================================
#ifdef DOI8B
  subroutine same_shape_pixels_nest_8(nside, templ8, list, reflexion, nrep)
    use num_rec, only : iindexx
    integer(i4b), parameter :: MKD = I8B
#else
  subroutine same_shape_pixels_nest  (nside, templ8, list, reflexion, nrep)
    use num_rec, only : iindexx
    integer(i4b), parameter :: MKD = I4B
#endif

    character(len=*), parameter :: code = "same_shape_pixels_nest"
    integer(I4B),    intent(IN)                        :: nside
    integer(MKD),    intent(IN)                        :: templ8
    integer(MKD),    pointer,   dimension(:), optional :: list
    integer(I4B),    pointer,   dimension(:), optional :: reflexion
    integer(MKD),    intent(OUT)            , optional :: nrep

    integer(MKD) :: ntplt, nr, i
    integer(MKD),   allocatable,   dimension(:)   :: idx, tmp
    integer(i4b) :: ll, lr
    logical(lgt) :: do_list, do_rfx, do_nrep
    !-------------------------------------------------------------

    do_list = (present(list))
    do_rfx  = (present(reflexion))
    do_nrep = (present(nrep))
  
    if (do_rfx .and. .not. do_list) then
       print*,'Error in '//code//'. Can not have Reflexion without pixel list'
       call fatal_error
    endif

    ntplt= nside2ntemplates(nside)

    if (templ8 < 0 .or. templ8 >= ntplt) then
       print*,'Error on template argument'
       print*,'Nside = ',nside,', Template = ',templ8
       print*,'Template should be in [0, (1+Nside*(Nside+6))/4-1=',ntplt-1,']'
       call fatal_error(code//' Abort')
    endif

    ! makes pixel list in RING indexing
    call same_shape_pixels_ring( nside, templ8, list, reflexion, nr)
    if (do_nrep) nrep = nr

    if (do_list .or. do_rfx) then
       allocate(idx  (0:nr-1))
       allocate(tmp  (0:nr-1))

       ! convert to NESTED
       if (do_list) then
          ll = lbound(list,      dim=1)
          do i = 0, nr-1
             call ring2nest(nside, list(ll+i), list(ll+i))
          enddo

          ! sort arrays in increasing order of pixel NESTED index
          call iindexx(nr, list, idx)
          !!!!!  list(ll:nr-1+ll) = list(idx-1+ll) ! replaced by loops for Gfortran
          do i=0, nr-1
             tmp(i) = list(idx(i)-1+ll)
          enddo
          do i=0, nr-1
             list(ll+i) = tmp(i)
          enddo
       endif
  
       if (do_rfx) then
          lr = lbound(reflexion,      dim=1)
          !!!!! reflexion(lr:nr-1+lr) = reflexion(idx-1+lr) ! replaced by loops for Gfortran
          do i=0, nr-1
             tmp(i) = reflexion(idx(i)-1+lr)
          enddo
          do i=0, nr-1
             reflexion(lr+i) = tmp(i)
          enddo
       endif

       deallocate(tmp)
       deallocate(idx)
    endif

    return

#ifdef DOI8B
  end subroutine same_shape_pixels_nest_8
#else
  end subroutine same_shape_pixels_nest
#endif


  !=======================================================================
  ! uniq2nest
  !  returns the Nside and Pnest (NESTED ordering pixel index) of
  !  input Puniq (Unique pixel indexing)
  !   Puniq = Pnest + 4 * Nside^2
  ! Sept 2015
  !=======================================================================

#ifdef DOI8B
  subroutine uniq2nest_8(puniq, nside, pnest)
    integer(i4b), parameter :: MKD = I8B
#else
  subroutine uniq2nest  (puniq, nside, pnest)
    integer(i4b), parameter :: MKD = I4B
#endif
    integer(MKD), intent(in)  :: puniq
    integer(I4B), intent(out) :: nside
    integer(MKD), intent(out) :: pnest
    integer(MKD) :: tmp_i
    real(DP)     :: tmp_d

    tmp_i = cheap_isqrt(puniq / 4_MKD)
    tmp_d = log( real(tmp_i, kind=DP) ) / log(2.0_DP)
    nside = 2_i4b ** int( tmp_d, kind=I4B )
    pnest = puniq - (4_MKD * nside) * nside

    return
#ifdef DOI8B
  end subroutine uniq2nest_8
#else
  end subroutine uniq2nest
#endif

  !=======================================================================
  ! nest2uniq
  !  returns Puniq (Unique pixel indexing)
  !  for given Nside and Pnest (NESTED ordering pixel index) of
  !   Puniq = Pnest + 4 * Nside^2
  ! Sept 2015
  !=======================================================================

#ifdef DOI8B
  subroutine nest2uniq_8(nside, pnest, puniq)
    integer(i4b), parameter :: MKD = I8B
#else
  subroutine nest2uniq  (nside, pnest, puniq)
    integer(i4b), parameter :: MKD = I4B
#endif
    character(len=*), parameter :: code='nest2uniq'
    integer(I4B), intent(in)  :: nside
    integer(MKD), intent(in)  :: pnest
    integer(MKD), intent(out) :: puniq
    integer(MKD) :: npix

    npix = nside2npix(nside)
    call assert(npix > 0, code//": invalid Nside")
    puniq = pnest + (4_MKD * nside) * nside


#ifdef DOI8B
    end subroutine nest2uniq_8
#else
    end subroutine nest2uniq
#endif

  !=======================================================================
  !=======================================================================
  !=======================================================================
! #ifdef DOI8B
!     integer(i4b), parameter :: MKD = I8B
! #else
!     integer(i4b), parameter :: MKD = I4B
! #endif
! #ifdef DOI8B
!     integer(i4b), parameter :: MKD = I8B
! #else
!     integer(i4b), parameter :: MKD = I4B
! #endif

! #ifdef DOI8B
! #else
! #endif
