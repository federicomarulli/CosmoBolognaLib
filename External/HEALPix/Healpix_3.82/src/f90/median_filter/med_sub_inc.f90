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
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Computes the computes the median filtered map from a full sky map in HEALPIX pixelisation
! Written and developed by E. Hivon (efh@ipac.caltech.edu) 
! based on HEALPIX pixelisation developed by K. Gorski
!
! Copyright 1997-2005 by Eric Hivon and Krzysztof M. Gorski.
!  All rights reserved.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!=======================================================================
  !=======================================================================
  !  EXTERNAL LIBRARY:
  !     this code uses the CFITSIO library that can be found at
  !     http://heasarc.gsfc.nasa.gov/docs/software/fitsio/
  !
  !  RELATED LITTERATURE:
  !     about HEALPIX   : see Gorski et al, 2005, ApJ, 622, 759
  !     about this code : see Hivon & Gorski, in preparation
  !
  !  HISTORY:
  !     April 2005, Eric Hivon, IPAC
  !
  !  FEEDBACK:
  !     for any questions : hivon@iap.fr
  !               Jun 2010 : supports large maps
  !
  !=======================================================================
  !     version 2.2
  !=======================================================================
  ! this file can not be compiled on its own.
  ! It must be inserted into the file median_filter.f90 by the command  include
  !
  integer(I4B),             parameter :: KCL   = KMAP  ! precision of c(l) arrays
  !-------------------------------------------------------------------------

  integer(I4B) :: nsmax

  real(kind=KMAP),     DIMENSION(:,:),   ALLOCATABLE :: map_in
  real(kind=KMAP),     DIMENSION(:,:),   ALLOCATABLE :: map_mf
  integer(kind=I4B) :: status

  integer(kind=I8B) :: npixtot
  integer(kind=I4B) :: nmaps, n_pols
  integer(kind=I4B) :: ordering, order_type
  integer(kind=I4B) :: i
  integer(kind=I4B) :: simul_type, polar_fits, polar
  integer(kind=i4b) :: type

  real(kind=KMAP),    parameter :: fbad_value = HPX_DBADVAL
  real(kind=DP)                 :: pix_size_arcmin
  real(kind=DP)                 :: mf_radius_arcmin, mf_radius_rad, mf_radius_deg
  real(kind=SP)                 :: clock_time, time0, time1
  real(kind=SP)                 :: ptime, ptime0, ptime1

  character(len=80), DIMENSION(1:120) :: header
  integer(kind=I4B) :: nlheader

  character(len=FILENAMELEN)          :: infile, mffile
  character(len=FILENAMELEN)          :: healpixtestdir
  character(len=FILENAMELEN)          :: description
  character(len=100)                  :: chline, sstr
  LOGICAL(kind=LGT) :: polarisation, fill_holes
  character(len=*), parameter :: VERSION = HEALPIX_VERSION

  type(paramfile_handle) :: handle

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
  write(*,'(a)') " *** Median filter for a Temperature/Polarisation map ***"
  if (KMAP == SP) print*,'                Single precision outputs'
  if (KMAP == DP) print*,'                Double precision outputs'
  PRINT *, " "
  handle = parse_init(parafile)

  !     --- choose temp. only or  temp. + pol. ---
  description = concatnl( &
       & " Do you want to filter", &
       & " 1) a Temperature only map", &
       & " 2) Temperature + polarisation maps")
  simul_type = parse_int(handle, 'simul_type', default=1, vmin=1, vmax=2, descr=description)
  if (simul_type .eq. 1) polarisation = .false.
  if (simul_type .eq. 2) polarisation = .true.

  !     --- gets the file name for the map ---
  chline = "''"
  healpixtestdir = get_healpix_test_dir()
  if (trim(healpixtestdir)/='') chline = trim(healpixtestdir)//'/map.fits'
  description = concatnl( &
       & "", &
       & " Enter input file name (Map FITS file): ")
  infile = parse_string(handle, 'infile',default=chline, descr=description, filestatus='old')

  !     --- finds out the pixel number of the map and its ordering ---
  npixtot = getsize_fits(infile, nmaps = nmaps, ordering=ordering, nside=nsmax, type = type, polarisation = polar_fits)
  if (nsmax<=0) then
     print*,"Keyword NSIDE not found in FITS header!"
     call fatal_error(code)
  endif
  if (type == 3) npixtot = nside2npix(nsmax) ! cut sky input data set
  if (nsmax/=npix2nside(npixtot)) then
     print*,"FITS header keyword NSIDE does not correspond"
     print*,"to the size of the map!"
     call fatal_error(code)
  endif

  if (polarisation .and. (nmaps >=3) .and. polar_fits == -1) then
     print*,"The input fits file MAY NOT contain polarisation data."
     print*,"Proceed at your own risk"
  endif

  if (polarisation .and. (nmaps<3 .or. polar_fits ==0)) then
     print *,"The file does NOT contain polarisation maps"
     print *,"only the temperature field will be filtered"
     polarisation = .false.
  endif

  !     --- check ordering scheme ---
  if ((ordering/=1).and.(ordering/=2)) then
     PRINT*,"The ordering scheme of the map must be RING or NESTED."
     PRINT*,"No ordering specification is given in the FITS-header!"
     call fatal_error(code)
  endif
  
  order_type = ordering

  ! ----------  default parameters ----------------
  npixtot = nside2npix(nsmax)
  pix_size_arcmin = 360.*60./SQRT(npixtot*PI)
  !   -------------------------------------
  write(sstr,'(g20.5)') pix_size_arcmin

  !     --- gets the disk radius  ---
  description = concatnl(&
       & "", &
       & " Pixel size in arcmin: "//trim(sstr), &
       & " Enter the radius in ARCMIN of the disk on which the median is to be computed")
  mf_radius_arcmin = parse_double(handle, 'mf_radius_arcmin', &
       &                       vmin=0.0_dp, default=3*pix_size_arcmin, descr=description)
  mf_radius_deg = mf_radius_arcmin / 60.0_dp
  mf_radius_rad = mf_radius_deg * DEG2RAD


  description = concatnl(&
       & "", &
       & " How do you want to treat missing pixels:",&
       & " - replace them by the median of surrounding disc (.true./t/yes/y/1) ", &
       & " - leave them unchanged               (fill_holes=.false./f/no/n/0)")
  fill_holes = parse_lgt(handle, 'fill_holes', default=.false., descr=description)
  
  !     --- gets the output median filtered filename ---
  description = concatnl(&
       & "", &
       & " Enter output file name for filtered map (eg: !medfilt.fits): ")
  mffile = parse_string(handle, 'mffile', default= '!medfilt.fits', descr=description, filestatus='new')

  PRINT *," "
  call parse_check_unused(handle, code=lcode)
  call parse_summarize(handle,code=lcode,prec=KMAP)
  call parse_finish(handle)
  call brag_openmp()

  !-----------------------------------------------------------------------
  !              allocate space for arrays
  !-----------------------------------------------------------------------
  polar = 0
  if (polarisation) polar = 1
  n_pols = 1 + 2*polar ! either 1 or 3

  ALLOCATE(map_in(0:npixtot-1,1:n_pols),stat = status)
  call assert_alloc(status,code,"map_in")

  ALLOCATE(map_mf(0:npixtot-1,1:n_pols),stat = status)
  call assert_alloc(status,code,"map_mf")

  !-----------------------------------------------------------------------
  !                      reads the map
  !-----------------------------------------------------------------------
  PRINT *,"      "//code//"> Inputting map "
  call input_map(infile, map_in, npixtot, n_pols, &
       & fmissval=fbad_value, header=header, ignore_polcconv=.true.)

  !-----------------------------------------------------------------------
  !                      filter the map
  !-----------------------------------------------------------------------
  PRINT *,"      "//code//"> Filtering map "
  do i=1, n_pols
     call medfiltmap(map_in(0:,i), mf_radius_rad, map_mf(0:,i), &
          &          nest=order_type-1_i4b, fmissval=fbad_value, fill_holes=fill_holes)
  enddo
  !-----------------------------------------------------------------------
  !                      output the map
  !-----------------------------------------------------------------------
  PRINT *,"      "//code//"> Outputting map "

  ! update header
  call add_card(header,"COMMENT", "-----------------------------------------------")
  call add_card(header,"HISTORY","median filtered map")
  call add_card(header,"HISTORY","input map: "//trim(infile))
  call add_card(header,"CREATOR",code,    update = .true.)
  call add_card(header,"VERSION",version, update = .true.)
  call add_card(header,"MFRADIUS",mf_radius_deg,"[Deg] median filter radius")
  call add_card(header,"BAD_DATA",fbad_value,"Sentinel value given to missing data")
  call add_card(header,"COMMENT","-----------------------------------------------")
  
  ! call output_map(map_mf, header, mffile) !not on Sun
  nlheader = size(header)
  call write_bintab(map_mf, npixtot, n_pols, header, nlheader, mffile)

  deallocate( map_in )
  deallocate( map_mf )

  !-----------------------------------------------------------------------
  !                      report card
  !-----------------------------------------------------------------------
  call wall_clock_time(time1)
  call cpu_time(ptime1)
  clock_time = time1 - time0
  ptime      = ptime1 - ptime0

  write(*,9000) " "
  write(*,9000) " Report Card for "//code//" run"
  write(*,9000) "----------------------------------------"
  write(*,9000) " "
  if (.not. polarisation) then
     write(*,9000) "      Temperature alone"
  else
     write(*,9000) "    Temperature + Polarisation"
  endif
  write(*,9000) " "
  write(*,9000) " Input map              : "//TRIM(infile)
  write(*,9010) " Number of pixels       : ", npixtot
  write(*,9020) " Pixel size in arcmin   : ", pix_size_arcmin
  write(*,9020) " Filter radius in arcmin: ", mf_radius_arcmin
  write(*,9000) " Filtered map           : "//TRIM(mffile)
  write(*,9030) " Clock and CPU time [s] : ", clock_time, ptime

  !-----------------------------------------------------------------------
  !                       end of routine
  !-----------------------------------------------------------------------

  write(*,9000) " "
  write(*,9000) " "//code//"> normal completion"

9000 format(a)
9005 format(a,i8)
9010 format(a,i16)
9020 format(a,g20.5)
9030 format(a,f11.2,f11.2)

