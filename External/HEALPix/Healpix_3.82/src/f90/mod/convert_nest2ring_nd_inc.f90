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
! generic body of the subroutines convert_nest2ring_*_nd 
! to be inserted as is in pix_tools.f90
!    edited 2009-04-03 to accept Nside>8192
!--------------------------------------------------------------
    !=======================================================================
    !     makes the conversion NEST to RING
    !=======================================================================
    integer(kind=I4B), intent(in) :: nside
    integer(kind=I4B) :: nd, j
    integer(kind=I4B) :: ipn4, ipr4
    integer(kind=I8B) :: ipn8, ipr8, npix
    integer(kind=I4B),  dimension(:), allocatable :: mapping4
    integer(kind=I8B),  dimension(:), allocatable :: mapping8
    !=======================================================================

    npix = nside2npix(nside)
    call assert (npix>0,       code//": invalid Nside")

    nd = size(map,2)

    ! case where 2nd dimension=1
    if (nd == 1) then
       map1 => map(0:npix-1,1)
       call convert_nest2ring(nside, map1)
       return
    endif

    allocate(map_tmp(0:npix-1))

    if (nside <= ns_max4) then
       allocate(mapping4(0:npix-1))

    ! do N->R mapping only once
!$OMP parallel default(none) &
!$OMP   shared(mapping4, npix, nside) private(ipr4, ipn4)
!$OMP do schedule(dynamic,64)
       do ipn4 = 0_i4b, npix-1
          call nest2ring(nside, ipn4, ipr4)
          mapping4(ipn4) = ipr4
       enddo
!$OMP end do
!$OMP end parallel

       ! convert maps one by one
       do j = 1, nd
          do ipn4 = 0_i4b, npix-1
             map_tmp(mapping4(ipn4)) = map(ipn4,j)
          enddo
          map(0:npix-1, j) = map_tmp(0:npix-1)
       enddo

       deallocate(mapping4)

    else

       allocate(mapping8(0:npix-1))

    ! do N->R mapping only once
!$OMP parallel default(none) &
!$OMP   shared(mapping8, npix, nside) private(ipr8, ipn8)
!$OMP do schedule(dynamic,64)
       do ipn8 = 0_i8b, npix-1
          call nest2ring(nside, ipn8, ipr8)
          mapping8(ipn8) = ipr8
       enddo
!$OMP end do
!$OMP end parallel

       ! convert maps one by one
       do j = 1, nd
          do ipn8 = 0_i8b, npix-1
             map_tmp(mapping8(ipn8)) = map(ipn8,j)
          enddo
          map(0:npix-1, j) = map_tmp(0:npix-1)
       enddo

       deallocate(mapping8)

    endif

    deallocate(map_tmp)

    return
