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
! generic body of the subroutines convert_inplace_*_nd 
! to be inserted as is in pix_tools.f90
!
! 2009-06-16: adapted for Nside > 8192
!--------------------------------------------------------------
    integer(kind=I1B), parameter :: TODO = 0_I1B, DONE = 1_I1B
    !
    external  subcall ! required by some compilers (DEC, Portland, ...)
    integer(kind=i4b)          :: nside, nd, j
    integer(kind=i8b)          :: npix
    integer(kind=i8b)          :: npix_8, ilast_8, i1_8, i2_8
    integer(kind=i4b)          :: npix_4, ilast_4, i1_4, i2_4
    integer(kind=I1B), dimension(:), allocatable::check
    !------------------------------------------------------------------
    npix=long_size(map,1)
    nd  =long_size(map,2)
    nside=npix2nside(npix)
    call assert (nside<=NS_MAX, code//": map too big")
    call assert (nside>0,       code//": invalid Nside")
    call assert (nd<=ND_MAX,    code//": map 2nd dim too large")

    print*, "Convert: Converting map pixelisation scheme"
    allocate(check(0:npix-1))
    check = TODO

    if (nside <= ns_max4) then 

       npix_4 = npix
       ilast_4=0                   !start at first pixel
       do
          pixbuf2(1:nd) = map(ilast_4,1:nd)      !initialise
          i1_4 = ilast_4
          call subcall(nside, i1_4, i2_4)
          do
             if (check(i2_4) == DONE) exit
             do j = 1, nd ! should be faster than 2 separate (implicit) loops
                pixbuf1(j) = map(i2_4,j)
                map(i2_4,j)  = pixbuf2(j)
             enddo
             check(i2_4) = DONE
             pixbuf2(1:nd) = pixbuf1(1:nd)
             i1_4 = i2_4
             call subcall(nside, i1_4, i2_4)
          enddo
          do
             if (.not. (check(ilast_4)==DONE .and. ilast_4<npix_4-1)) exit ! npix-1 or npix
             ilast_4 = ilast_4 + 1
          enddo
          if(ilast_4 == npix_4-1) exit ! npix-1 or npix
       enddo

    else

       npix_8 = npix
       ilast_8=0                   !start at first pixel
       do
          pixbuf2(1:nd) = map(ilast_8,1:nd)      !initialise
          i1_8 = ilast_8
          call subcall(nside, i1_8, i2_8)
          do
             if (check(i2_8) == DONE) exit
             do j = 1, nd ! should be faster than 2 separate (implicit) loops
                pixbuf1(j) = map(i2_8,j)
                map(i2_8,j)  = pixbuf2(j)
             enddo
             check(i2_8) = DONE
             pixbuf2(1:nd) = pixbuf1(1:nd)
             i1_8 = i2_8
             call subcall(nside, i1_8, i2_8)
          enddo
          do
             if (.not. (check(ilast_8)==DONE .and. ilast_8<npix_8-1)) exit ! npix-1 or npix
             ilast_8 = ilast_8 + 1
          enddo
          if(ilast_8 == npix_8-1) exit ! npix-1 or npix
       enddo

    endif

    deallocate(check)
    return

