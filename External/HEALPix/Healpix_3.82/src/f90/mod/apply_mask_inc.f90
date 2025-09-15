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
!--------------------------------------------------------------
!
! generic body of the subroutines apply_mask_*
! to be inserted as is in pix_tools.F90
!
! 2018-10-22: creation
! 2018-11-09: same behavior as remove_dipole, select_rings 
!             and libsharp's hpsharp_make_healpix_geom_info_2
!
!--------------------------------------------------------------

  real(KMAP), dimension(0:,1:), intent(inout) :: map
  integer(I4B),                 intent(in)    :: order
  real(KMAP), dimension(0:,1:), intent(in), optional    :: mask
  real(DP),   dimension(1:2),   intent(in), optional    :: zbounds

  real(DP) :: theta1, theta2, frac
  integer(I8B) :: npix_map, n_maps, i, p, k
  integer(I8B) :: npix_mask, n_masks, ncpix, ncp8
  integer(I8B),   dimension(:), allocatable :: cut_pixel8
  integer(I4B) :: nside
  character(len=*), parameter :: code = 'apply_mask'
  logical(LGT) :: verbose = .false.
!=====================================================================

  npix_map  = long_size(map,  1)
  n_maps    = long_size(map,  2)


  ! apply mask
  if (present(mask)) then
     npix_mask = long_size(mask, 1)
     n_masks   = long_size(mask, 2)
     if (npix_map == npix_mask) then 
        do i=1,n_maps
           k = min(i, n_masks) ! if single mask, apply it to TQU
           do p=0,npix_map-1
              map(p,i) = map(p,i) * mask(p,k)
           enddo
        enddo
     else
        if (verbose) then
           print*,npix_map, npix_mask
           print*,code//': mask is not applied'
        endif
     endif
  endif

  ! apply zbounds
  if (present(zbounds)) then
     !!!!if (abs(zbounds(1)+1.0_dp) + abs(zbounds(2)-1.0_dp) > 1.d-15) then ! not (-1,1)
     if ((zbounds(1)+1.0_dp)> EPSILON(1.0_DP) .or. (zbounds(2)-1.0_dp) < -EPSILON(1.0_DP)) then ! not (-1,1)
        nside = npix2nside(npix_map)
        if (zbounds(2) <= zbounds(1)) then ! remove one strip
           frac = (zbounds(1)-zbounds(2)) * 0.5_dp
        else ! remove 2 strips
           frac = (2.0_dp + zbounds(1)-zbounds(2)) * 0.5_dp
        endif
        ncpix  = min(nint(npix_map * frac * 1.1_dp + 10*nside,I8B), npix_map)
        theta1 = acos(zbounds(1)) 
        theta2 = acos(zbounds(2))
        allocate(cut_pixel8(0:ncpix))
        call query_strip(nside, theta1, theta2, cut_pixel8, ncp8, nest=order-1, inclusive=0)
        do i=1,n_maps
           map(cut_pixel8(0:ncp8-1),i) = 0.0_KMAP
        enddo
        deallocate(cut_pixel8)
     else
        if (verbose) then
           print*,code//': zbounds is ignored'
        endif
     endif
  endif

  return

