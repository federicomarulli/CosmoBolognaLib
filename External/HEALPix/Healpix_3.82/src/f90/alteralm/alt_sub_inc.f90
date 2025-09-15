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
!
! facility to modify a set of a_lm (spherical harmonics coefficients)
!
! Copyright 1997 by Eric Hivon and Krzysztof M. Gorski.
!  All rights reserved.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!=======================================================================
  !=======================================================================
  !  EXTERNAL LIBRARY:
  !     this code uses the FITSIO library that can be found at
  !     http://heasarc.gsfc.nasa.gov/docs/software/fitsio/fitsio.html
  !
  !  RELATED LITTERATURE:
  !     about HEALPIX   : see Gorski et al, 2004
  !     about this code : see Hivon & Gorski, 1997, in preparation
  !
  !  HISTORY:
  !     Dec 2004, Eric Hivon, IPAC
  !
  !  FEEDBACK:
  !     for any questions : hivon@iap.fr
  !
  !=======================================================================
  !     version 2.2
  !=======================================================================
  ! this file can not be compiled on its own.
  ! It must be inserted into the file smoothing.f90 by the command  include
  !
  integer(I4B),             parameter :: KALM  = KMAP ! precision of alm arrays
  integer(I4B),             parameter :: KALMC = KIND((1.0_KMAP, 1.0_KMAP)) ! precision of alm arrays
  !-------------------------------------------------------------------------
  integer(I4B) :: nsmax_in, nlmax_in, nmmax_in, nalms
  integer(I4B) :: nsmax_out, nlmax_out, nmmax_out, nalms_out, nlmin_out
  integer(I4B) :: npol, ncl, nlheader, junk
  integer(I4B) :: i, k, l, m, status, cnt, lmm

  complex(KALMC), dimension(:,:,:), allocatable, target :: alm_sig, alm_err
  complex(KALMC), dimension(:,:),   pointer             :: palms
  real   (KALM),  dimension(:,:,:), allocatable         :: alms_in, alms_out
  real   (DP),    dimension(:,:),   allocatable         :: pixel_in, pixel_out
  real   (DP),    dimension(:,:),   allocatable         :: beam_in, beam_out
  real   (DP),    dimension(:,:),   allocatable         :: window_out

  real   (KALM)                       :: fwhm_deg
  real     (DP)                       :: fwhm_arcmin_in, fwhm_arcmin_out, fwhm_arcmin_def
  character(len=FILENAMELEN)          :: infile_alms,    outfile_alms
  character(len=FILENAMELEN)          :: beam_file_in,   beam_file_out, beam_file_def
  character(len=80), dimension(:,:), allocatable :: header
 !  character(len=*),         parameter :: CODE = 'ALTERALM'
  character(len=*), parameter :: VERSION = HEALPIX_VERSION
  character(len=10)                   :: sstr
  real     (SP) :: clock_time, time0, time1
  real     (SP) :: ptime, ptime0, ptime1

  character(len=20) :: coord_in, coord_out
  character(len=1), allocatable, dimension(:) :: coord_options
  real(dp) :: phi1, phi2, theta
  real(dp) :: epoch_in, epoch_out

  logical(LGT) :: do_alm_err, do_rot

  integer(I4B)      :: ktype, nsmax
  character(len=20) :: suffix, tag
  character(len=FILENAMELEN) :: def_dir, def_file, usr_dir, usr_file, final_file
  character(len=FILENAMELEN) :: windowname, windowfile, windowfile_in, windowfile_out
  logical(LGT) :: ok

  type(paramfile_handle)              :: handle
!  character(len=FILENAMELEN)          :: parafile = ''
  character(len=FILENAMELEN)          :: description
  !======================================================================

  call wall_clock_time(time0)
  call cpu_time(ptime0)
  !     --- read parameters interactively if no command-line arguments
  !     --- are given, otherwise interpret first command-line argument
  !     --- as parameter file name and read parameters from it:

  !     --- announces program, default banner ***
  print *, " "
  print *,"                    "//CODE//" "//version
  write(*,'(a)') &
       & " *** Modification of a set of alm (spherical harmonics coefficients) ***"
  if (KMAP == SP) print*,'                Single precision outputs'
  if (KMAP == DP) print*,'                Double precision outputs'
  print *, " "
  handle = parse_init(parafile)

  description = concatnl(&
       & "", &
       & " Enter the FITS file containing the input a_lm coefficients")
  infile_alms = parse_string(handle, 'infile_alms', descr=description,filestatus='old')
  call read_par(infile_alms, nsmax_in, nlmax_in, ncl, mmax=nmmax_in)
  if (nmmax_in < 0) nmmax_in = nlmax_in
  nalms = number_of_alms(infile_alms, extnum = npol)
  junk = getsize_fits(infile_alms, fwhm_arcmin = fwhm_arcmin_def, beam_leg = beam_file_def,&
       &             coordsys = coord_in, extno=0)
  do_alm_err = (ncl == 5)

  !     *********** FWHM or beam file ************
  !           ---------   IN   ---------
  if (fwhm_arcmin_def < 0) then
     description = concatnl(&
          & "", &
          & " The beam FWHM could not be found in the input a_lm file, ", &
          & " Enter the Gaussian beam FWHM to correct the a_lm from, ", &
          & "  in arcmin  >= 0 (eg: 0.) ")
     fwhm_arcmin_def = 0.0_dp
  else
     description = concatnl(&
          & "", &
          & " Enter the Gaussian beam FWHM to correct the a_lm from, ", &
          & "  in arcmin  >= 0 (eg: 0.) ")
  endif
  fwhm_arcmin_in = parse_double(handle, 'fwhm_arcmin_in', default=fwhm_arcmin_def, vmin=0.0_dp,descr=description)

  beam_file_in = ''
  description = concatnl(&
       & "", &
       & " Enter an external file name containing ", &
       & " the window function to correct the alm from (eg: mybeam.fits)", &
       & " NB: if set to an existing file, it will override the FWHM chosen above", &
       & "     if set to '', the gaussian FWHM will be used instead.")
  beam_file_in = parse_string(handle, 'beam_file_in', default=beam_file_def, &
       &                 descr=description, filestatus='old')
  if (trim(beam_file_in) /= '') then
     fwhm_arcmin_in = 0.
     print*,'fwhm_arcmin_in is now : 0.'
     print*,'The beam file '//trim(beam_file_in)//' will be used instead.'
  endif

  if (fwhm_arcmin_in < 0.0 .and. trim(beam_file_in) == '') then
     call fatal_error(CODE//"> no information found on input alms beam")
  endif

  !           ---------   OUT   ---------
  description = concatnl(&
       & "", &
       & " Enter the Gaussian beam FWHM to apply to the alm, ", &
       & "  in arcmin  >= 0 (eg: 0.) ")
  fwhm_arcmin_out = parse_double(handle, 'fwhm_arcmin_out', default=fwhm_arcmin_in, vmin=0.0_dp,descr=description)

  beam_file_out = ''
  description = concatnl(&
       & "", &
       & " Enter an external file name containing ", &
       & " the window function to apply to the alm (eg: mybeam.fits)", &
       & " NB: if set to an existing file, it will override the FWHM chosen above", &
       & "     if set to '', the gaussian FWHM will be used.")
  beam_file_out = parse_string(handle, 'beam_file_out', default="''", &
       &                 descr=description, filestatus='old')
  if (trim(beam_file_out) /= '') then
     fwhm_arcmin_out = 0.
     print*,'fwhm_arcmin_out is now : 0.'
     print*,'The beam file '//trim(beam_file_out)//' will be used instead.'
  endif
  fwhm_deg = fwhm_arcmin_out / 60.0_dp

  !     *********** coordinate system ************
  !           ---------   IN   ---------
  epoch_in = 2000.0_dp
  epoch_out = 2000.0_dp
  allocate(coord_options(1:9))
  coord_options = (/ 'C', 'E', 'Q', 'G', 'c', 'e', 'q', 'g', ' ' /)
  coord_in=''
  if (.not. any(coord_options(1:8)==trim(adjustl(coord_in)))) then
     description = concatnl(&
       & "", &
       & " The astronomical coordinate system of the input file was not found, ", &
       & "  or is not valid: "//trim(coord_in), &
       & " If you want to ROTATE the alm, you will have to provide one as below,", &
       & " if not, leave it blank", &
       & " G = Galactic, E = Ecliptic, C/Q = Celestial = eQuatorial")
     coord_in = parse_string(handle,'coord_in',default='',descr=description,options=coord_options)
  else
     description = concatnl(&
       & "", &
       & " If you want to ROTATE the alm, provide the input astro. coordinates as below,", &
       & " if not, leave it blank", &
       & " G = Galactic, E = Ecliptic, C/Q = Celestial = eQuatorial")
     coord_in = parse_string(handle, 'coord_in',default=coord_in,descr=description, options=coord_options)
  endif
  description = concatnl(&
       & "", &
       & " Provide the coordinate epoch used for Ecliptic precession")
  epoch_in = parse_double(handle, 'epoch_in',default=2000.0_dp,descr=description)

  !           ---------   OUT   ---------
  coord_out=''
  if (trim(coord_in) /= '' .or. handle%interactive) then
     description = concatnl(&
          & "", &
          & " Provide the output astronomical coordinate system")
     coord_out = parse_string(handle, 'coord_out',default=coord_in,descr=description, options=coord_options)
  endif
  description = concatnl(&
       & "", &
       & " Provide the coordinate epoch used for Ecliptic precession")
  epoch_out = parse_double(handle, 'epoch_out',default=epoch_in,descr=description)
  do_rot = (trim(coord_in) /= '' .and. trim(coord_in) /= trim(coord_out))
  do_rot = (do_rot .or. (epoch_in /= epoch_out))

  !     *********** pixel size ************
  !           ---------   IN   ---------
  if (nsmax_in < 0) then
     description = concatnl(&
          & "", &
          & " The resolution parameter (Nside) of the input file was not found ", &
          & " Enter the correct one (or 0 if the alm were free of pixel window)")
     nsmax_in = parse_int(handle, 'nsmax_in', vmin=0, descr=description)
  endif

  !           ---------   OUT   ---------
  sstr = trim(string(nsmax_in,'(i6)'))
  description = concatnl(&
       & "", &
       & " The resolution parameter (Nside) for the alm was "//sstr, &
       & " Enter the Nside whose pixel window function will be applied to the a_lm, ",&
       & " after correcting for the input one ",&
       & " (Nside HAS to be a power of 2, eg: 32, 512, ..., or 0 for no pixel)" )
  nsmax_out = parse_int(handle, 'nsmax_out', default=nsmax_in, vmin=0, descr=description)


  !     --- check for pixel-window-files ---
  do ktype=1,2
     if (ktype == 1) then
        nsmax = nsmax_in
        suffix = '_in'
        tag = ' (INPUT alm) '
     endif
     if (ktype == 2) then
        nsmax = nsmax_out
        suffix = '_out'
        tag = ' (OUTPUT alm) '
     endif

     if (nsmax > 0) then
        windowname = get_healpix_pixel_window_file(nsmax)
        def_file = trim(windowname)
        def_dir  = get_healpix_data_dir()

22      continue
        final_file = ''
        ok = .false.
        ! if interactive, try default name in default directories first
        if (handle%interactive) ok = scan_directories(def_dir, def_file, final_file)
        if (.not. ok) then
           ! not found, ask the user
           description = concatnl("",&
                &        " Could not find "//trim(tag)//" window file", &
                &        " Enter the directory where this file can be found:")
           usr_dir = parse_string(handle,'winfiledir'//trim(suffix),default="''",descr=description)
           if (trim(usr_dir) == '') usr_dir = trim(def_dir)
           description = concatnl("",&
                &        " Enter the name of the "//trim(tag)//" window file:")
           usr_file = parse_string(handle,'windowfile'//trim(suffix),default=def_file,descr=description)
           ! look for new name in user provided or default directories
           ok   = scan_directories(usr_dir, usr_file, final_file)
           ! if still fails, crash or ask again if interactive
           if (.not. ok) then
              print*,' File not found'
              if (handle%interactive) goto 22
              call fatal_error(code)
           endif
        endif
        windowfile = final_file
     else
        windowfile = ''
     endif
     if (ktype == 1) windowfile_in  = windowfile
     if (ktype == 2) windowfile_out = windowfile
  enddo

  !     *********** multipole range ************
  !           ---------   IN   ---------
  nlmin_out = 0
  sstr = trim(string(nlmax_in,'(i6)'))
!   lmm = min(3*nsmax_out-1, nlmax_in)
  lmm = nlmax_in
  description = concatnl(&
       & "", &
       & " The maximum multipole for the alm was "//sstr, &
       & " Enter the maximum multipole for the output." )
  nlmax_out = parse_int(handle, 'nlmax_out', default=nlmax_in, vmin=0, vmax=lmm, descr=description)
  nmmax_out = min(nlmax_out, nmmax_in)

  !           ---------   OUT   ---------
  description = concatnl(&
       & "", &
       & " Enter the FITS file to contain the output a_lm coefficients")
  outfile_alms = parse_string(handle, 'outfile_alms', descr=description,filestatus='new')

  print *," "
!  call parse_summarize(handle)
  call parse_check_unused(handle, code=lcode)
  call parse_summarize(handle,code=lcode,prec=KMAP)
  call parse_finish(handle)
  call brag_openmp()

  !------------------------------------
  ! allocate arrays
  !------------------------------------
  nlheader = 120
  allocate(header(1:nlheader, 1:npol))
  allocate(alm_sig(1:npol, 0:nlmax_in, 0:nmmax_in), stat = status)
  call assert_alloc(status,CODE,"alms_sig")
  if (do_alm_err) then
     allocate(alm_err(1:npol, 0:nlmax_in, 0:nmmax_in), stat = status)
     call assert_alloc(status,CODE,"alms_err")
  endif

  !------------------------------------
  ! input alm's
  !------------------------------------
  allocate(alms_in(0:nalms-1,1:ncl+1,1:npol),stat = status)
  call assert_alloc(status,CODE,"alms_in")

  print *,"      "//CODE//"> Inputting alms "
  !!  print*,trim(infile_alms),nalms,ncl,nlheader,npol
  call fits2alms(infile_alms, nalms, alms_in, ncl, header, nlheader, npol)
  !!  print*,'alms read'
  do k = 1, npol ! loop on T, [E, B]
     do i = 0, nalms-1
        l = nint(alms_in(i,1,k))
        m = nint(alms_in(i,2,k))
        if (m < 0) then
           write(*,'(a,i6,a,i6)') '     the a_lm with m<0 are ignored, l=',l,', m=',m
        else if (l > nlmax_in .or. m> nmmax_in .or. m>l) then
           print*,'Error on l,m',l,m
           print*,'larger than ',nlmax_in, nmmax_in
           call fatal_error(CODE)
        else
           alm_sig(k,l,m) = CMPLX(alms_in(i,3,k), alms_in(i,4,k), kind=KALM)
           if (do_alm_err) then 
              alm_err(k,l,m) = CMPLX(alms_in(i,5,k), alms_in(i,6,k), kind=KALM)
           endif
        endif
     enddo
  enddo
  deallocate(alms_in)

  !------------------------------------
  ! alter alm's
  !------------------------------------
  print *,"      "//CODE//"> Altering alms "

  ! generate input window function (pixel & beam)
  allocate(pixel_in(0:nlmax_in, 1:npol), beam_in(0:nlmax_in, 1:npol), stat=status)
  call assert_alloc(status, CODE, "pixel_in, beam_in")
  if (nsmax_in == 0) then
     call pixel_window(pixel_in, nside=nsmax_in)
  else
     call pixel_window(pixel_in, windowfile=windowfile_in)
  endif
  call generate_beam(fwhm_arcmin_in, nlmax_in, beam_in, beam_file=beam_file_in)

  ! generate output window function (pixel & beam)
  allocate(pixel_out(0:nlmax_out, 1:npol), beam_out(0:nlmax_out, 1:npol), stat=status)
  call assert_alloc(status, CODE, "pixel_out, beam_out")
  if (nsmax_out == 0) then
     call pixel_window(pixel_out, nside=nsmax_out)
  else
     call pixel_window(pixel_out, windowfile=windowfile_out)
  endif
  call generate_beam(fwhm_arcmin_out, nlmax_out, beam_out, beam_file=beam_file_out)

  ! generate effective window function
  allocate(window_out(0:nlmax_out, 1:npol), stat=status)
  call assert_alloc(status, CODE, "window_out")

  window_out = pixel_in(0:nlmax_out, 1:npol) * beam_in(0:nlmax_out, 1:npol)
  where(window_out == 0.0) window_out = 1.e-15
  window_out = (pixel_out * beam_out) / window_out

  call alter_alm(0_i4b, nlmax_out, nmmax_out, 0.0_KALM, alm_sig, &
       &            window=real(window_out,kind=KALM))
  if (do_alm_err) then
     call alter_alm(0_i4b, nlmax_out, nmmax_out, 0.0_KALM, alm_err, &
          &            window=real(window_out,kind=KALM))
  endif

  ! rotate alm's
  if (do_rot) then
     call brag_openmp()
     coord_in  = trim(strupcase(coord_in))
     coord_out = trim(strupcase(coord_out))
9700 format(7x,a,a," (",f6.1,")",a,a," (",f6.1,")")
     write(*,9700) CODE//'> Rotating alms:   ', &
          & trim(coord_in),epoch_in,' -> ',trim(coord_out),epoch_out
     ! figure out Euler angles corresponding to coordinate changes
     call coordsys2euler_zyz(epoch_in, epoch_out, coord_in, coord_out, phi1, theta, phi2)
!     print*,phi1,theta,phi2
     ! perform the rotation
     call rotate_alm(nlmax_out, alm_sig, phi1, theta, phi2)
     if (do_alm_err) then
        call rotate_alm(nlmax_out, alm_err, phi1, theta, phi2)
     endif
  endif
  !------------------------------------
  ! output alm's
  !------------------------------------
  print *,"      "//CODE//"> Outputting alms "

  ! update header read from input file
  do k = 1, npol
     if (trim(beam_file_out) == "") then
        call add_card(header(1:,k),"FWHM"    ,fwhm_deg   ," [deg] FWHM of gaussian symmetric beam")
     else
        call add_card(header(1:,k),"BEAM_LEG",trim(beam_file_out),"File containing Legendre transform of symmetric beam")
     endif
     call add_card(header(1:,k),"NSIDE", nsmax_out, "Resolution parameter for HEALPIX")
     call add_card(header(1:,k),"MAX-LPOL",nlmax_out  ,"Maximum multipole l ") ! updated by FITS writting routines
     call add_card(header(1:,k),"MIN-LPOL",nlmin_out  ,"Minimum multipole l ") ! updated by FITS writting routines
     call add_card(header(1:,k),"MAX-MPOL",nmmax_out  ,"Maximum m ")           ! updated by FITS writting routines
     call add_card(header(1:,k),"COORDSYS",coord_out  ,"Coordinates system ")
     call add_card(header(1:,k),"EPOCH",   epoch_out  ,"Coordinates system epoch")
     call add_card(header(1:,k),"HISTORY"," file "//trim(infile_alms))
     call add_card(header(1:,k),"HISTORY"," modified by "//CODE//" "//VERSION)
     if (fwhm_arcmin_out /= fwhm_arcmin_in) call add_card(header(1:,k),"HISTORY"," changed FWHM")
     if (beam_file_out /= beam_file_in) call add_card(header(1:,k),"HISTORY"," changed BEAM_LEG")
     if (nsmax_out /= nsmax_in) call add_card(header(1:,k),"HISTORY"," changed NSIDE")
     if (nlmax_out /= nlmax_in) call add_card(header(1:,k),"HISTORY"," changed MAX-LPOL")
     if (coord_out /= coord_in) call add_card(header(1:,k),"HISTORY"," changed COORDSYS")
     if (epoch_out /= epoch_in) call add_card(header(1:,k),"HISTORY"," changed EPOCH")
  enddo

  if (do_alm_err) then
     nalms_out = ((nmmax_out+1)*(2*nlmax_out+2-nmmax_out))/2
     allocate(alms_out(0:nalms_out-1, 1:ncl+1, 1:npol), stat=status)
     call assert_alloc(status, CODE, "alms_out")
     do k = 1, npol
        cnt = 0
        do l = 0, nlmax_out
           do m = 0, min(l, nmmax_out)
              alms_out(cnt, 1, k) = l
              alms_out(cnt, 2, k) = m
              alms_out(cnt, 3, k) = real (alm_sig(k,l,m), kind=KALM) !       ---> col 2
              alms_out(cnt, 4, k) = aimag(alm_sig(k,l,m))            !       ---> col 3
              alms_out(cnt, 5, k) = real (alm_err(k,l,m), kind=KALM) !       ---> col 4
              alms_out(cnt, 6, k) = aimag(alm_err(k,l,m))            !       ---> col 5
              cnt = cnt + 1
           enddo
        enddo
        call assert(cnt <= nalms_out, CODE//' Error on index')
     enddo
     call alms2fits(outfile_alms, nalms_out, alms_out, ncl, header, nlheader, npol)
     deallocate(alms_out)
  else
     do k = 1, npol
        palms => alm_sig(k,0:nlmax_out,0:nmmax_out)
        call dump_alms (outfile_alms, palms, nlmax_out, header(:,k), nlheader, k-1_i4b)
     enddo
  endif
  call wall_clock_time(time1)
  call cpu_time(ptime1)
  clock_time = time1 - time0
  ptime      = ptime1 - ptime0

  write(*,9000) " "
  write(*,9000) " Report Card for "//code//" run"
  write(*,9000) "------------------------------------"
  write(*,9000) " "
  write(*,9000) " Input a_lms           : "//TRIM(infile_alms)
  write(*,9005) " Input Nside           : ",nsmax_in
  write(*,9005) " Input multipole range : 0 < l <= ", nlmax_in
  write(*,9300) " Input coordinates     : ",trim(coord_in),epoch_in
  if (trim(beam_file_in) == '') then
     write(*,9020) " Input Gauss. FWHM ['] : ", fwhm_arcmin_in
  else
     write(*,9000) " Input beam file       : "//trim(beam_file_in)
  endif
  write(*,9000) " Output a_lms          : "//TRIM(outfile_alms)
  if (nsmax_in /= nsmax_out) &
       & write(*,9005)       " Output Nside          : ",nsmax_out
  if (nlmax_in /= nlmax_out) &
       & write(*,9005) " Output multipole range: 0 < l <= ", nlmax_out
  if (do_rot) write(*,9300) " Output coordinates    : ",trim(coord_out),epoch_out
  if (trim(beam_file_out) == '') then
     write(*,9020) " Output Gauss. FWHM [']: ", fwhm_arcmin_out
  else
     write(*,9000) " Output beam file      : "//trim(beam_file_out)
  endif
  write(*,9030) " Clock and CPU time [s]: ", clock_time, ptime

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
9300 format(a,a," (",f6.1,")")
