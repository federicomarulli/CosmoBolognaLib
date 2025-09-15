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
!! Upgrade / Degrade a Healpix map
!! Written and developed by E. Hivon (efh@ipac.caltech.edu) and K. Gorski
!! (krzysztof.m.gorski@jpl.nasa.gov) based on HEALPIX pixelisation developed
!! by K. Gorski
!!
!! Copyright 1997 by Eric Hivon and Krzysztof M. Gorski.
!!  All rights reserved.
!!
!! author Eric Hivon (hivon@iap.fr)
!!
!  RELATED LITTERATURE:
!     about HEALPIX   : see Gorski et al, 2005, ApJ, 622, 759
!-----------------------------------------------------------------------
  !=======================================================================
  !     version 2.2
  !=======================================================================
  ! this file can not be compiled on its own.
  ! It must be inserted into the file ud_grade.f90 by the command  include
  !
  INTEGER(I4B) :: nside_in, nside_out

  REAL(KMAP),     DIMENSION(:,:),   ALLOCATABLE :: map_in
  REAL(KMAP),     DIMENSION(:,:),   ALLOCATABLE :: map_out
  INTEGER(I4B)  status

  INTEGER(I8B) npix_in, npix_out
  INTEGER(I4B) nmaps, type
  INTEGER(I4B) j, polar_fits

  CHARACTER(LEN=filenamelen)          ::  outfile, infile
  CHARACTER(LEN=80), DIMENSION(1:120) :: header, header_in
  CHARACTER(LEN=80) :: char
  INTEGER(I4B) :: nlheader
  INTEGER(I4B) :: ordering, order_type
!  REAL(KMAP)   :: fmissval = real(HPX_DBADVAL,kind=KMAP)
  REAL(KMAP)   :: fmissval = HPX_DBADVAL
  real(SP)                 :: clock_time, time0, time1
  real(SP)                 :: ptime, ptime0, ptime1
!  real(SP)                 :: ct1,ct2,ct,pt1,pt2,pt

  integer(i4b), parameter :: nm_max = 5
  CHARACTER(len=80), dimension(1:nm_max) :: ttype, com_ttype, tunit, com_tunit
  INTEGER(I4B), dimension(1:nm_max) :: ct_ttype, ct_tunit
  character(len=1), dimension(1:nm_max) :: sn = (/ '1', '2', '3', '4', '5' /)
  character(len=20) :: coordsys

!   CHARACTER(LEN=*), PARAMETER :: code = 'UD_GRADE'
  character(len=*), parameter :: VERSION = HEALPIX_VERSION

  type(paramfile_handle) :: handle
  character(len=20) :: polcconv
  character(len=80) :: com_polcconv
  integer(i4b)      :: n_polcconv

  !-----------------------------------------------------------------------
  !                    get input parameters and create arrays
  !-----------------------------------------------------------------------

  call wall_clock_time(time0)
  call cpu_time(ptime0)
  !     --- read parameters interactively if no command-line arguments
  !     --- are given, otherwise interpret first command-line argument
  !     --- as parameter file name and read parameters from it:

  !     --- announces program, default banner ***
  PRINT *, " "
  PRINT *, "                "//code//" "//version
  write(*,'(a)') '          *** Upgrading / Degrading of a map ***    '
  if (KMAP == SP) print*,'                Single precision outputs'
  if (KMAP == DP) print*,'                Double precision outputs'
  PRINT *, " "
  handle = parse_init(parafile)

  !     --- gets the file name for the map ---
  infile=parse_string(handle,'infile',default='map.fits', &
    descr='Input file name (Map FITS file, eg: map_CDM.fits):', &
    filestatus='old')
  PRINT *," "

  !     --- finds out the pixel number of the map and its ordering ---

  npix_in = getsize_fits(infile, nmaps=nmaps, ordering=ordering, nside=nside_in, &
       type=type, polarisation = polar_fits, coordsys = coordsys)
  if (nside_in.eq.0) then
     print*,'Keyword NSIDE not found in FITS header!'
     stop 1
  endif
  if (type /= 2) then  ! accept only full sky map
     print*,'> '//code//version//' only accepts full sky maps'
     print*,'> '//trim(infile)//' is not a full sky map'
     stop 1
  endif
  if (nside_in.ne.npix2nside(npix_in)) then ! accept only full sky map
     print*,'FITS header keyword NSIDE = ',nside_in,' does not correspond'
     print*,'to the size of the map!'
     stop 1
  endif


  !     --- gets the resolution of the output map ---
  WRITE(*,'(a,i5)') ' The input map has Nside = ',nside_in
  nside_out=parse_int(handle,'nside_out',default=64_i4b, &
    descr='Resolution parameter (Nside) for the output skymap:', &
    vmin=1)

  npix_out = nside2npix(nside_out)

  !     --- check ordering scheme ---

  if ((ordering.ne.1).and.(ordering.ne.2)) then
     PRINT*,'The ordering scheme of the map must be RING or NESTED.'
     PRINT*,'No ordering specification is given in the FITS-header!'
     stop 1
  endif

  order_type = ordering

  !     --- gets the output sky map filename ---
  outfile=parse_string(handle,'outfile',default='outmap.fits', &
    descr='Output map file name (eg. map_up.fits):', &
    filestatus='new')
  PRINT *," "
  call parse_check_unused(handle, code=lcode)
  call parse_summarize(handle,code=lcode,prec=KMAP)
  call parse_finish(handle)
  call brag_openmp()

  !-----------------------------------------------------------------------
  !              allocate space for arrays
  !-----------------------------------------------------------------------

  ALLOCATE(map_in(0:npix_in-1,1:nmaps),stat = status)
  call assert_alloc(status,code,'map_in')

  ALLOCATE(map_out(0:npix_out-1,1:nmaps),stat = status)
  call assert_alloc(status,code,'map_out')

  !-----------------------------------------------------------------------
  !                      reads the map
  !-----------------------------------------------------------------------
  PRINT *,'      '//code//'> Input original map '
  call input_map(infile, map_in(0:,1:nmaps), npix_in, nmaps, &
       &   fmissval=fmissval, header=header_in, ignore_polcconv=.true.)
  ! get input type and unit
  do j=1,nmaps
     call get_card(header_in, 'TTYPE'//sn(j), ttype(j), com_ttype(j), count=ct_ttype(j))
     call get_card(header_in, 'TUNIT'//sn(j), tunit(j), com_tunit(j), count=ct_tunit(j))
  enddo

  !-----------------------------------------------------------------------
  !                 makes the UP/DE-GRADATION
  !-----------------------------------------------------------------------
  if (nside_out .gt. nside_in) then
     PRINT *,'      '//code//'> UPgrade the map '
  else
     PRINT *,'      '//code//'> DEgrade the map '
  endif

!   call wall_clock_time(pt1)
!   call cpu_time(ct1)
  if (order_type == 1) then
     call udgrade_ring(map_in , nside_in, map_out, nside_out)
  else
     call udgrade_nest(map_in , nside_in, map_out, nside_out)
  endif
!   call wall_clock_time(pt2)
!   call cpu_time(ct2)
!   print*,pt2-pt1,ct2-ct1

  !-----------------------------------------------------------------------
  !                        generates header
  !-----------------------------------------------------------------------
  PRINT *,'      '//code//'> Writing up/de-graded map to FITS file '
  call write_minimal_header(header, 'map', &
       order = order_type, nside = nside_out, coordsys = coordsys, &
       creator = code, version = version, polar = (polar_fits == 1))
  if (nside_out .gt. nside_in) then
     call add_card(header,'EXTNAME','''UPGRADED DATA''', update = .true.)
  else
     call add_card(header,'EXTNAME','''DEGRADED DATA''', update = .true.)
  endif
  call add_card(header,'HISTORY','Input Map in '//TRIM(infile))
  write(char,'(i6)') nside_in
  call add_card(header,'HISTORY','Input Map resolution NSIDE =  '//TRIM(char))

  ! replicate in output file the POLCCONV found in input file (added 2017-02-27)
  call get_card(header_in, 'POLCCONV', polcconv, com_polcconv, count=n_polcconv)
  if (n_polcconv == 0) then ! delete what write_minimal_header wrote
     call del_card(header,'POLCCONV')
  else ! overwrite what write_minimal_header wrote
     call add_card(header,'POLCCONV', polcconv, update=.true.)
  endif

  ! copy input type and units
  do j=1,nmaps
     if (ct_ttype(j) > 0) then
        call add_card(header,'TTYPE'//sn(j), ttype(j), com_ttype(j), update = .true.)
     else
        call add_card(header,'TTYPE'//sn(j), 'unknown','unidentified type', update = .true.)
     endif
     if (ct_tunit(j) > 0) then
        call add_card(header,'TUNIT'//sn(j), tunit(j), com_tunit(j), update = .true.)
     else
        call add_card(header,'TUNIT'//sn(j), 'unknown','unknown units', update = .true.)
     endif
  enddo
  call add_card(header) ! blank line
  nlheader = SIZE(header)
  !-----------------------------------------------------------------------
  !                      write the map to FITS file
  !-----------------------------------------------------------------------
  call write_bintab(map_out, npix_out, nmaps, header, nlheader, outfile)

  !-----------------------------------------------------------------------
  !                      deallocate memory for arrays
  !-----------------------------------------------------------------------
  DEALLOCATE(map_in)
  DEALLOCATE(map_out)

  !-----------------------------------------------------------------------
  !                      output and report card
  !-----------------------------------------------------------------------
  call wall_clock_time(time1)
  call cpu_time(ptime1)
  clock_time = time1 - time0
  ptime      = ptime1 - ptime0

  WRITE(*,9000) " "
  WRITE(*,9000) " Report Card for "//code//" run"
  WRITE(*,9000) " -----------------------------"
  WRITE(*,9000) " "

  WRITE(*,9000) " Input map              : "//TRIM(infile)
  WRITE(*,9005) " Number of maps         : ", nmaps
  WRITE(*,9010) " Number of pixels (IN)  : ", npix_in
  WRITE(*,9000) " Output map             : "//TRIM(outfile)
  WRITE(*,9010) " Number of pixels (OUT) : ", npix_out
  write(*,9030) " Clock and CPU time [s] : ", clock_time, ptime
  !-----------------------------------------------------------------------
  !                       end of routine
  !-----------------------------------------------------------------------
  WRITE(*,9000) " "
  WRITE(*,9000) " "//code//"> normal completion"

9000 format(a)
9005 format(a,i8)
9010 format(a,i16)
9020 format(a,g20.5)
9030 format(a,f11.2,f11.2)

