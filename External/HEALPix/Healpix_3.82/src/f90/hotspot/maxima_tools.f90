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
module maxima_tools
   ! This module contains
   !    find_maxima     : performs optimised search of maxima and minima
   !                      inside a face and uses neighbours_nest
   !                      to find the neighbours of a pixel on the border
   !                      of the face.
   ! Benjamin D. Wandelt October 1997
   ! 2010-11-23, E. Hivon: accept maps with Nside > 8192
   !
   use healpix_types
   use pix_tools, only : npix2nside, xy2pix_nest, neighbours_nest
   implicit none
contains
   !====================================================================
   subroutine find_maxima(map,face_num,peak)
      !====================================================================
      !     Finds all maxima in the face face_num
      !     INPUTS:
      !     map:      REAL array, length NPIX which contains the map
      !     face_num: INTEGER Number of the face in which extremum search is to be
      !               performed
      !     OUTPUTS:  
      !     peak:     INTEGER array, length NPIX which contains 1 in the pixels
      !               where maxima were found and -1 in the pixels where minima
      !               were found.
      !
      !     Benjamin D. Wandelt October 1997
      ! 2010-11-23, E. Hivon: accept maps with Nside > 8192
      !====================================================================
      use bit_manipulation
      implicit none

      real(sp), intent(in) :: map(0:)
      integer(i4b), intent(in):: face_num
      integer(i4b), intent(out)::peak(0:)

      integer(i8b) :: npix, nsidesq, origin
      integer(i8b) :: ipix0m,ipix00,ipix0p,ipixpm,ipixp0,ipixpp
      integer(i8b) :: neighb(8)

      integer(i4b) :: nside, nn
      integer(i4b) :: ix,iy,ixp,iyp,iym
      real(sp) :: pixval0m,pixval00,pixval0p, &
           &     pixvalpm,pixvalp0,pixvalpp, &
           &     pixvalmm,pixvalm0,pixvalmp
      real(sp) :: comp
!      real(sp), intrinsic :: sqrt

      !-----------------------------------------------------------------
      if (face_num<0.or.face_num>11) &
           &     stop 'find_maxima: face_num out of range'
      npix=size(map)
      ! if(npix>12*ns_max*ns_max) stop "find_maxima: map too big"
      nside   = npix2nside(npix)
      nsidesq = nside * nside

      !     First do interior
      !     Reduce number of calls to xy2pix_nest
      !     by reusing values from previous step. This brings 9 calls
      !     down to 3. What about parallelisation of inner loop?

      do iy=0+1,nside-1-1
         ix=0
         ixp=ix+1
         iyp=iy+1
         iym=iy-1
         call xy2pix_nest(nside, ix , iym, face_num, ipix0m)
         call xy2pix_nest(nside, ix , iy , face_num, ipix00)
         call xy2pix_nest(nside, ix , iyp, face_num, ipix0p)
         call xy2pix_nest(nside, ixp, iym, face_num, ipixpm)
         call xy2pix_nest(nside, ixp, iy , face_num, ipixp0)
         call xy2pix_nest(nside, ixp, iyp, face_num, ipixpp)
         pixval0m = map(ipix0m)
         pixval00 = map(ipix00)
         pixval0p = map(ipix0p)
         pixvalpm = map(ipixpm)
         pixvalp0 = map(ipixp0)
         pixvalpp = map(ipixpp)
         do ix=0+1,nside-1-1
            ipix00=ipixp0
            pixvalmm=pixval0m
            pixvalm0=pixval00
            pixvalmp=pixval0p
            pixval0m=pixvalpm
            pixval00=pixvalp0
            pixval0p=pixvalpp
            ixp=ix+1
            call xy2pix_nest(nside, ixp, iym, face_num, ipixpm)
            call xy2pix_nest(nside, ixp, iy , face_num, ipixp0)
            call xy2pix_nest(nside, ixp, iyp, face_num, ipixpp)
            pixvalpm = map(ipixpm)
            pixvalp0 = map(ipixp0)
            pixvalpp = map(ipixpp)
            comp=max(pixvalmm,pixvalm0,pixvalmp, &
                 &           pixval0m,pixval0p,pixvalpm,pixvalp0,pixvalpp)
            if(pixval00.gt.comp)then
               peak(ipix00)=+1 ! Maximum
            endif
            comp=min(pixvalmm,pixvalm0,pixvalmp, &
                 &           pixval0m,pixval0p,pixvalpm,pixvalp0,pixvalpp)
            if(pixval00.lt.comp)then
               peak(ipix00)=-1 ! Minimum
            endif
         enddo
      enddo

      !----------------------
      !     Now do boundaries
      !----------------------

      origin=face_num*nsidesq
      iy=0
      !*$*assert concurrent call
      !c$doacross local(ipix00,pixval00,neighb,nn,comp),
      !c$& share(origin,face_num,peak,map,nside)
      do ix=0, nside-1
         call xy2pix_nest(nside,ix,iy,face_num, ipix00)
         pixval00=map(ipix00)
         call neighbours_nest(nside,ipix00,neighb,nn)
         comp=maxval(map(neighb(1:nn)))
         if(pixval00.gt.comp)then
            peak(ipix00)=+1 !Maximum
         endif
         comp=minval(map(neighb(1:nn)))
         if(pixval00.lt.comp)then
            peak(ipix00)=-1 !Minimum
         endif
         if(ix.gt.0)then
            ipix00=origin+mod(swapLSBMSB(ipix00),nsidesq) !E-W flip
            pixval00=map(ipix00)
            call neighbours_nest(nside,ipix00,neighb,nn)
            comp=maxval(map(neighb(1:nn)))
            if(pixval00.gt.comp)then
               peak(ipix00)=+1 !Maximum
            endif
            comp=minval(map(neighb(1:nn)))
            if(pixval00.lt.comp)then
               peak(ipix00)=-1 !Minimum
            endif
         endif
      enddo

      ix=nside-1
      !*$*assert concurrent call
      !c$doacross
      !c local(ipix00,pixval00,neighb,nn,comp), share(peak)
      do iy=1, nside-1
         call xy2pix_nest(nside,ix,iy,face_num, ipix00)
         pixval00=map(ipix00)
         call neighbours_nest(nside,ipix00,neighb,nn)
         comp=maxval(map(neighb(1:nn)))
         if(pixval00.gt.comp)then
            peak(ipix00)=+1 !Maximum
         endif
         comp=minval(map(neighb(1:nn)))
         if(pixval00.lt.comp)then
            peak(ipix00)=-1 !Minimum
         endif
         if(iy.lt.(nside-1))then
            ipix00=origin+mod(swapLSBMSB(ipix00),nsidesq) !E-W flip
            pixval00=map(ipix00)
            call neighbours_nest(nside,ipix00,neighb,nn)
            comp=maxval(map(neighb(1:nn)))
            if(pixval00.gt.comp)then
               peak(ipix00)=+1 !Maximum
            endif
            comp=minval(map(neighb(1:nn)))
            if(pixval00.lt.comp)then
               peak(ipix00)=-1 !Minimum
            endif
         endif
      enddo

      return
   end subroutine find_maxima

end module maxima_tools
