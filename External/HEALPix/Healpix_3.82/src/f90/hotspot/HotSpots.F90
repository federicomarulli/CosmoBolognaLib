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
!======================================================================
program HotSpots
  !======================================================================
  !     Finds extrema in the 2-D field on the sphere.
  !     Uses file "param.dat" for inputs
  !
  ! Benjamin D. Wandelt October 1997
  ! 2010-11-23, E. Hivon: accept maps with Nside > 8192
  !======================================================================
  use healpix_types
  use maxima_tools, ONLY : find_maxima
  use fitstools,  ONLY : getsize_fits, input_map, write_bintab
#ifdef G95
  use pix_tools,  ONLY : convert_ring2nest, convert_nest2ring, npix2nside
#else
  use pix_tools,  ONLY : convert_inplace, ring2nest, nest2ring, npix2nside
#endif
  USE head_fits,  ONLY : add_card, write_minimal_header
  use misc_utils, ONLY : assert_alloc, fatal_error
  USE extension,  ONLY : getArgument, nArguments
  use paramfile_io, ONLY : paramfile_handle, parse_init, parse_string
  use long_intrinsic, only: long_count

  implicit none
  integer(i4b) :: npix, nmin, nmax, i
  integer(i4b) :: nside
  integer(i4b) :: face_num

  real(sp), pointer ::  map(:)
  real(sp), allocatable :: inmap(:,:)
  real(sp), allocatable :: realpeak(:,:)
  integer(i4b),pointer :: peak(:)

  integer(i4b) :: ordering,order_type
  integer(i4b) :: nlheader
  integer(i4b) :: nmaps = 1

  character(len=80), dimension(1:120) :: header
  character(len=filenamelen) :: parafile,infile,outfile_extrema, &
                                outfile_max,outfile_min
  character(len=20) :: coordsys
  integer(i4b) :: status

  type(paramfile_handle) :: handle

  CHARACTER(LEN=*), PARAMETER :: code    = 'HOTSPOTS'
  character(len=*), parameter :: VERSION = HEALPIX_VERSION

  !--------------------------------------------------------------------------
100 FORMAT(A)

  !-----------------------------------------------------------------------
  !                    get input parameters and create arrays
  !-----------------------------------------------------------------------

  !     --- read parameters interactively if no command-line arguments
  !     --- are given, otherwise interpret first command-line argument
  !     --- as parameter file name and read parameters from it:
  if (nArguments() == 0) then
     parafile=''
  else
     if (nArguments() /= 1) then
        print '(" Usage: hotspots [parameter file name]")'
        call fatal_error(code)
     end if
     call getArgument(1,parafile)
  end if

  handle=parse_init (parafile)

  !     --- announces program, default banner ***
  PRINT *,'                    '//code//' '//version
  PRINT *,' ***  Finding Maxima and Minima of a HEALPIX Map  ***'
  PRINT *,' *** (Application of the Neighbour-finding Tools) ***'
  PRINT *,''

  !     --- gets the file name for the map ---
  infile=parse_string(handle,'infile',default='map.fits', &
       descr='Input file name (Map FITS file, e.g. map_CDM.fits):', &
       filestatus='old')
  PRINT *,''

  !     --- finds out the pixel number of the map and its ordering ---
  npix = getsize_fits(infile ,ordering=ordering, nside=nside, coordsys = coordsys)

  if (nside.eq.0) then
     print*,'Keyword NSIDE not found in FITS header!'
     call fatal_error(code)
  endif
  if (nside.ne.npix2nside(npix)) then
     print*,'FITS header keyword NSIDE does not correspond'
     print*,'to the size of the map!'
     call fatal_error(code)
  endif

  if(nside<=1) then
     PRINT*,' '//code//'> Nside needs to be at least 2 for extrema calculations.'
     call fatal_error(code)
  end if

  !     --- check ordering scheme ---

  if ((ordering.ne.1).and.(ordering.ne.2)) then
     PRINT*,'The ordering scheme of the map must be RING or NESTED.'
     PRINT*,'No ordering specification is given in the FITS-header!'
     call fatal_error()
  endif

  order_type = ordering

!     --- gets the output extrema map filename ---
  outfile_extrema=parse_string(handle,'extrema_outfile', &
       default='pixlminmax.fits', &
       descr='Extrema map file name (eg, pixlminmax.fits) :', &
       filestatus='new')
  PRINT *,''

!     --- gets the output filename for ascii file with MAXIMA ---
  outfile_max=parse_string(handle,'maxima_outfile', &
       default='maxima.dat', &
       descr='Maxima file name (eg, maxima.dat) :')
!       filestatus='new')
  PRINT *,''

!     --- gets the output filename for ascii file with MINIMA ---
  outfile_min=parse_string(handle,'minima_outfile', &
       default='minima.dat', &
       descr='Maxima file name (eg, minima.dat) :')
!       filestatus='new')
  PRINT *,''

  allocate(map(0:npix-1),stat=status)
  call assert_alloc(status,code,'map')

  allocate(inmap(0:npix-1,1),stat=status)
  call assert_alloc(status,code,'inmap')

  PRINT*, ' '//code//'> Reading '//TRIM(infile)//'...'

  call input_map(infile,inmap,npix,nmaps)

  map(0:npix-1) = inmap(0:npix-1,1)
  deallocate(inmap)

  if(order_type==1) then
     write(*,*) code//'> Converting to NESTED numbering scheme.'
#ifdef G95
     call convert_ring2nest(nside, map)
#else
     call convert_inplace(ring2nest,map)
#endif
  endif

  allocate(peak(0:npix-1))
  call assert_alloc(status,code,'peak')
  peak=0

  write(*,*) code//'> Starting extrema search in pixel space.'

  ! This is fully parallelisable. All 12 faces can be searched synchronously
  !*$*assert concurrent call
  ! peak = +1 for maximum, -1 for minimum, 0 otherwise
  do face_num=0,11
     call find_maxima(map,face_num,peak)
  end do

  write(*,*) 'Done.'

  if(order_type==1) then
#ifdef G95
     write(*,*) code//'> Converting map back to RING numbering.'
     call convert_nest2ring(nside,map)
     write(*,*) code//'> Converting peaks to RING numbering.'
     call convert_nest2ring(nside,peak)
#else
     write(*,*) code//'> Converting map back to RING numbering.'
     call convert_inplace(nest2ring,map)
     write(*,*) code//'> Converting peaks to RING numbering.'
     call convert_inplace(nest2ring,peak)
#endif
  endif

  write (*,*) code//' Writing '//trim(outfile_extrema)//' FITS file...'
  call write_minimal_header(header, 'map', &
       nside = nside, order = ordering, coordsys = coordsys, &
       polar = .false., creator = code, version = version)

  call add_card(header,'EXTNAME','PEAK MAP', update = .true.)
  call add_card(header,'COMMENT', 'Input Map ='//TRIM(infile))
  call add_card(header,'COMMENT', 'Local Min and Max keep their original value, the rest is set to 0')
  call add_card(header,'TTYPE1', 'PEAKS','Maxima and Minima', update = .true.)

  ! ------
  ! outputs a FITS map set to 0 but for extrema that keep their original value
  allocate(realpeak(0:npix-1,1))
  realpeak=0
  where(peak/=0) realpeak(:,1)=map

  nlheader = SIZE(header)
  call write_bintab(realpeak, npix, 1_i4b, header, nlheader, outfile_extrema)

  deallocate(realpeak)

  ! outputs the set of maxima (pixel number and value)
  ! and the set of minima (pixel number and value) in 2 ascii files
  open(14,file=outfile_max,status="unknown")
  open(15,file=outfile_min,status="unknown")
  write(*,*) 'Output: Writing '//trim(outfile_max)//' and '//trim(outfile_min)//'.'

  nmax=long_count(peak>0)
  write(*,*) "Found ",nmax," Maxima"
  nmin=long_count(peak<0)
  write(*,*) "Found ",nmin," Minima"

  do i=0,npix-1
     if(peak(i)>0) then     !Maximum
        write(14,'(i10,f12.6)') i,map(i)
     elseif(peak(i)<0) then !Minimum
        write(15,'(i10,f12.6)') i,map(i)
     endif
  enddo

  close(14)
  close(15)

  deallocate(peak)
  deallocate(map)
  stop
end program HotSpots
