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
! Simulation of full sky map of CMB anisotropies from their power spectrum
! in HEALPIX pixelisation.
! Written and developed by E. Hivon (efh@ipac.caltech.edu) and K. Gorski
! (krzysztof.m.gorski@jpl.nasa.gov) based on HEALPIX pixelisation developed
! by K. Gorski
!
! Copyright 1997 by Eric Hivon and Krzysztof M. Gorski.
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
  !     about this code : see Hivon & Gorski, 1997, in preparation
  !
  !  HISTORY:
  !     April-October 1997, Eric Hivon, TAC
  !               Jan 1998 : translation in Fortran 90
  !               Aug 1998 : version 1.0.0
  !               Apr 1999 : constr.real, precomp.plms included, FKH
  !               Nov 2000 : v 1.2, EH
  !               Sep 2002 : implement new parser
  !               2004-2005: v 2.0
  !               2006-2007: v 2.1
  !               Jun 2010 : supports large maps
  !
  !  FEEDBACK:
  !     for any questions : hivon@iap.fr
  !
  !=======================================================================
  !     version 2.2.0
  !=======================================================================
  ! this file can not be compiled on its own.
  ! It must be inserted into the file synfast.f90 by the command  include
  !
  integer(I4B),             parameter :: KALM  = KMAP ! precision of alm arrays
  integer(I4B),             parameter :: KALMC = KIND((1.0_KMAP, 1.0_KMAP)) ! precision of alm arrays
  integer(I4B),             parameter :: KCL   = KMAP  ! precision of c(l) arrays
  !-------------------------------------------------------------------------

  integer(I4B) :: nsmax, nlmax, nmmax

  COMPLEX(KALMC), DIMENSION(:,:,:), ALLOCATABLE :: alm_TGC
  REAL(KALM),     DIMENSION(:,:),   ALLOCATABLE :: alms_in
  REAL(KMAP),     DIMENSION(:,:),   ALLOCATABLE :: map_TQU
  REAL(DP),       DIMENSION(:,:),   ALLOCATABLE :: plm
  real(DP),       dimension(:,:),   allocatable :: beam_const, pixel_const
  integer(I4B) :: status

  integer(I4B) l, m         !, alms
  integer(I8B) npixtot
  integer(I4B) i, j, junk
  integer(I4B) iseed, ioriginseed
  integer(I4B) simul_type
  integer(I4B) polar, polar_cl
  integer(I4B) nlheader, nalms, extnum_alms
  integer(I8B) :: n_plm
  integer(I4B) :: plm_nside, plm_lmax, plm_pol, plm_mmax

  REAL(KALM) :: fwhm_arcmin, fwhm_deg, pix_size_arcmin
  REAL(SP) :: clock_time, time0, time1
  REAL(SP) :: ptime, ptime0, ptime1
  REAL(DP) :: nullval

  character(len=FILENAMELEN)          :: infile
  character(len=FILENAMELEN)          :: outfile
!  character(len=FILENAMELEN)          :: parafile = ''
  character(len=FILENAMELEN)          :: infile_alms
  character(len=FILENAMELEN)          :: outfile_alms
  character(len=FILENAMELEN)          :: infile_plm
  character(len=FILENAMELEN)          :: windowfile
  character(len=FILENAMELEN)          :: windowname
  character(len=FILENAMELEN)          :: def_dir, def_file
  character(len=FILENAMELEN)          :: usr_dir, usr_file
  character(len=FILENAMELEN)          :: final_file
  character(len=FILENAMELEN)          :: healpixtestdir
  character(len=FILENAMELEN)          :: beam_file
  character(len=80), DIMENSION(1:180) :: header_PS
  character(len=80), DIMENSION(1:180) :: header
!  character(len=4)                    :: sstr
  LOGICAL(LGT) :: ok, bad, polarisation, input_cl, input_alms, do_map, output_alms, apply_windows
!   character(len=*), PARAMETER :: code = "SYNFAST"
  character(len=*), parameter :: VERSION = HEALPIX_VERSION
  character(len=80), allocatable, dimension(:) :: units_alm, units_map

  type(paramfile_handle) :: handle
  character(len=FILENAMELEN)          :: description
  character(len=100)                  :: chline
  type(planck_rng) :: rng_handle

  integer(i4b) :: deriv, n_pols, n_maps, n_clin
  character(len=20) :: coordsys = ' '
  real(DP), dimension(1:2) :: zbounds
!  real :: ttime0, ttime1, ttime2

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
  PRINT *,"                    "//code//" "//version
  write(*,'(a)') &
       & " *** Synthesis of a Temperature (and Polarisation) map from its power spectrum ***"
  if (KMAP == SP) print*,'                Single precision outputs'
  if (KMAP == DP) print*,'                Double precision outputs'
  PRINT *, " "
  handle = parse_init(parafile)

  !     --- choose temp. only or  temp. + pol. ---
  description = concatnl( &
       & " Do you want to simulate", &
       & " 1) Temperature only", &
       & " 2) Temperature + polarisation", &
       & " 3) Temperature + 1st derivatives", &
       & " 4) Temperature + 1st & 2nd derivatives", &
       & " 5) T+P + 1st derivatives", &
       & " 6) T+P + 1st & 2nd derivates")
  simul_type = parse_int(handle, 'simul_type', default=1, vmin=1, vmax=6, descr=description)
  deriv = 0 ! no derivatives
  polarisation = (simul_type == 2 .or. simul_type == 5 .or. simul_type == 6)
  if (simul_type == 3 .or. simul_type == 5) deriv = 1
  if (simul_type == 4 .or. simul_type == 6) deriv = 2

  !     --- gets the effective resolution of the sky map ---
3 continue
  description = concatnl( &
       & "", &
       & " Enter the resolution parameter (Nside) for the simulated skymap: ",&
       & " (Npix = 12*Nside**2, where Nside HAS to be a power of 2, eg: 32, 512, ...)" )
  nsmax = parse_int(handle, 'nsmax', default=32, descr=description)
  if (nside2npix(nsmax) < 0) then
     print *, " Error: nsmax is not a power of two."
     if (handle%interactive) goto 3
     call fatal_error(code)
  endif

  !     --- gets the L range for the simulation ---
  WRITE(chline,"(a,i5,a)") "We recommend: (0 <= l <= l_max <= ",3*nsmax-1,")"
  description = concatnl(&
       & "", &
       & " Enter the maximum l range (l_max) for the simulation. ", &
       & chline )
  nlmax = parse_int(handle, 'nlmax', default=2*nsmax, descr=description)

  ! ----------  default parameters ----------------
  nmmax   = nlmax
  npixtot = nside2npix(nsmax)
  pix_size_arcmin = 360.*60./SQRT(npixtot*PI)
  !   -------------------------------------

  !     --- gets the file name for the Power Spectrum ---
  chline = ''
  healpixtestdir = get_healpix_test_dir()
  if (trim(healpixtestdir)/='') chline = trim(healpixtestdir)//'/cl.fits'
  description = concatnl(&
       & "", &
       & " Enter input Power Spectrum filename", &
       & " If '' (empty single quotes) is typed, no file is read ", &
       & "    (then an alm_file should be entered later on)")
  infile = parse_string(handle, 'infile', default=chline, descr=description, filestatus='old')
  input_cl = (trim(infile) /= "")

  !     --- inputs the random seed ---
  description = concatnl(&
       & "", &
       & " Enter NON-ZERO random seed (eg: -1) : ", &
       & "  if set to zero, it will be over-riden by self generated seed. " )
  iseed = parse_int(handle, 'iseed', default=-1, descr=description)
!!  iseed = -ABS(iseed)

  !     --- inputs the FWHM or beam file ---
  beam_file = ''
  description = concatnl(&
       & "", &
       & " Enter the Gaussian beam FWHM in arcmin >= 0 (eg: 420.) : ")
  fwhm_arcmin = parse_double(handle, 'fwhm_arcmin', default=420.0_dp, vmin=0.0_dp,descr=description)

  description = concatnl(&
       & "", &
       & " Enter an external file name containing ", &
       & " a symmetric beam Legendre transform (eg: mybeam.fits)", &
       & " NB: if set to an existing file, it will override the FWHM chosen above", &
       & "     if set to '', the gaussian FWHM will be used.")
  beam_file = parse_string(handle, 'beam_file', default="''", &
       &                 descr=description, filestatus='old')
  if (trim(beam_file) /= '') then
     fwhm_arcmin = 0.
     print*,'fwhm_arcmin is now : 0.'
     print*,'The beam file '//trim(beam_file)//' will be used instead.'
  endif
  fwhm_deg = fwhm_arcmin/60.

  !     --- check for pixel-window-files ---
  windowname = get_healpix_pixel_window_file(nsmax)

  def_file = trim(windowname)
  def_dir  = get_healpix_data_dir()

22 continue
  final_file = ''
  ok = .false.
  ! if interactive, try default name in default directories first
  if (handle%interactive) ok = scan_directories(def_dir, def_file, final_file)
  if (.not. ok) then
     ! not found, ask the user
     description = concatnl("",&
          &        " Could not find window file", &
          &        " Enter the directory where this file can be found:")
     usr_dir = parse_string(handle,'winfiledir',default="''",descr=description)
     if (trim(usr_dir) == '') usr_dir = trim(def_dir)
     description = concatnl("",&
          &        " Enter the name of the window file:")
     usr_file = parse_string(handle,'windowfile',default=def_file,descr=description)
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

  !     --- constrained realisations ---
  description = concatnl( &
       & "", &
       & " Enter name of file with constraining a_lms :", &
       & "(eg, alms.fits;  if '' is typed, no file is read)", &
       & " These should EITHER already have the same beam and pixel as the final map,", &
       & "OR be totally free of beam and pixel window.", &
       & " See the alteralm facility if neither is true.")
  infile_alms = parse_string(handle, 'almsfile', default="''", &
       &                     descr=description, filestatus="old")
  input_alms = (trim(infile_alms) /= "")

  if (input_alms) then

     nalms=number_of_alms(infile_alms,extnum=extnum_alms)
     junk = getsize_fits(infile_alms, coordsys=coordsys)

  !      --- apply windows to constraining alms ---
     description = concatnl( &
          & "", &
          & " Do you want to use the constraining a_lms as they are, (F) ", &
          & "or apply first the beam and pixel window functions chosen previously ? (T)", &
          & " The latter assumes the a_lm to be free of any window functions,", &
          & " see the alteralm facility if these a_lm already have different window functions.")
     apply_windows = parse_lgt(handle, 'apply_windows', default=.false., descr=description)

  endif

  ! --- check that some form of input has been given ---
  if (.not. input_cl .and. .not. input_alms) then
     print *," "
     print *," "//code//"> ERROR: No input spectrum were given"
     print *," "//code//"> ERROR:   and the a_lm were not constrained"
     print *," "//code//"> ERROR: No map can be generated"
     call fatal_error("Abort execution")
  endif

  !          --- precomputed plms  ---
40 continue
  n_plm   = 0
  plm_pol = 0
  description = concatnl(&
       & "", &
       & " Enter name of file with precomputed p_lms :", &
       & "(eg, plm.fits; if '' is typed, no file is read)" )
  infile_plm = parse_string(handle, 'plmfile', default="''", descr=description, filestatus='old')

  if (trim(infile_plm) /= "") then
     ! check that the plm file matches the current simulation (for nside,lmax,mmax)
     call read_par(infile_plm, plm_nside, plm_lmax, plm_pol, mmax=plm_mmax)

     if ((plm_pol /= 1).and.(.not.polarisation)) plm_pol=1

     if (plm_nside /= nsmax .or. plm_lmax /= nlmax .or. plm_mmax /= nmmax) then
        print *," "//code//"> Plm file does not match map to simulate in Nside, Lmax or Mmax ! "
        print*,'nside (plm file, map) = ',plm_nside,nsmax
        print*,'lmax                  = ',plm_lmax,nlmax
        print*,'mmax                  = ',plm_mmax,nmmax
        if (handle%interactive) goto 40
        call fatal_error(code)
     endif
     n_plm = (nmmax + 1_i8b)*(2*nlmax + 2_i8b - nmmax)*nsmax
     print *," "
  endif

  !     --- gets the output sky map filename ---
  description = concatnl(&
       & "", &
       & " Enter Output map file name (eg, test.fits) :", &
       & "  (or !test.fits to overwrite an existing file)" , &
       & " (If '', no map is created) ")
  outfile = parse_string(handle, "outfile", &
       & default="!test.fits", descr=description, filestatus="new")
  do_map = (trim(outfile) /= '')

  !     --- gets the output alm-filename ---
  description = concatnl(&
       & "", &
       & " Enter file name in which to write a_lm used to synthesize the map : ", &
       & " (eg alm.fits or !alms.fits to overwrite an existing file): ", &
       & " (If '', the alms are not written to a file) ")
  outfile_alms = parse_string(handle, 'outfile_alms', &
       & default="''", descr=description, filestatus='new')
  output_alms = (trim(outfile_alms) /= '')

  PRINT *," "
  call parse_check_unused(handle, code=lcode)
  call parse_summarize(handle,code=lcode,prec=KMAP)
  call parse_finish(handle)
  call brag_openmp()

  !     --- If necessary use current system time to set seed  ---
  IF (iseed .eq. 0) THEN
     CALL SYSTEM_CLOCK(COUNT = iseed)
     !-- Seed should initially be set to a large, odd integer --
     IF (MOD(iseed,2) .EQ. 0) iseed = iseed + 1
     PRINT *,"      "//code//"> Generating random number seed ", iseed
  END IF
  ioriginseed = iseed

  ! --- initialize random sequence ---
  call rand_init(rng_handle, iseed) ! takes up to 4 seeds simultaneously

  !-----------------------------------------------------------------------
  !                  allocates space for arrays
  !-----------------------------------------------------------------------
  polar = 0
  if (polarisation) polar = 1
  n_pols = 1 + 2*polar ! either 1 or 3
  n_maps = max(n_pols, 3*deriv*n_pols) ! either 1 (T), 3 (TQU or T+d1) or 6 (T+d1+d2) or 9 (T,Q,U+d1) or 18 (T,Q,U+d1+d2)

  ALLOCATE(units_alm(1:n_pols),units_map(1:n_maps),stat = status)
  call assert_alloc(status,code,"units_alm & units_map")
  units_alm(:)=''

  ALLOCATE(alm_TGC(1:n_pols, 0:nlmax, 0:nmmax),stat = status)
  call assert_alloc(status,code,"alm_TGC")

  if (do_map) then
     ALLOCATE(map_TQU(0:npixtot-1, 1:n_maps),stat = status)
     call assert_alloc(status,code,"map_TQU")
  endif


  !-----------------------------------------------------------------------
  !                    creates/reads (and write) the a_lm
  !-----------------------------------------------------------------------

  PRINT *,"      "//code//"> Generating alm coefficients "
  PRINT *,"      "//code//"> Pixel window file:"//trim(adjustl(windowfile))

  alm_TGC = 0. ! very important

  if (input_cl) then
     polar_cl =     polar  ! 0 or 1: input 1 or 4 field spectrum
     !polar_cl = 2 * polar  ! 0 or 2; input 1 or 6 field spectrum
     !junk = getsize_fits(infile, nmaps=n_clin)
     !if (n_clin == 6) polar_cl = 2*polar ! input 1 or 6 spectra
     call create_alm(nsmax, nlmax, nmmax, polar_cl, &
          & infile, rng_handle, fwhm_arcmin, &
          & alm_TGC, header_PS, windowfile, units=units_alm, beam_file=beam_file)
     units_map(1:n_pols) = units_alm(1:n_pols)
     call del_card(header_PS, (/ "TUNIT#","TTYPE#"/))
  else
     header_PS = ""
     units_map = "unknown"
  endif

  !      --- constrained alms ---
  ! they are assumed to produce a REAL map
  ! therefore the a_lm with m<0 are ignored
  ! and are assumed to be a_l-|m| = (-1)^m a_l|m|
  if (input_alms) then

     allocate(alms_in(0:nalms-1,1:6),stat = status)
     call assert_alloc(status,code,"alms_in")
     alms_in = 0.0_KALM

     if (apply_windows) then
        allocate(pixel_const(0:nlmax, 1:n_pols), beam_const(0:nlmax, 1:n_pols), stat=status)
        call assert_alloc(status, CODE, "pixel_const, beam_const")
        call pixel_window(pixel_const, windowfile=windowfile)
        call generate_beam(real(fwhm_arcmin,kind=dp), nlmax, beam_const, beam_file=beam_file)
        beam_const = pixel_const * beam_const
        deallocate(pixel_const)
     else
        allocate(beam_const(0:nlmax, 1:n_pols), stat=status)
        beam_const = 1.0_KALM
     endif

     PRINT *,"      "//code//"> Reading alms for constrained realisations "
     if (apply_windows) PRINT *,"      "//code//"> and applying same beam and pixel size on them"
     do j = 1, MIN(n_pols,extnum_alms)
        call read_conbintab(infile_alms, alms_in(0:nalms-1,1:6), nalms,extno=j-1_i4b)
        do i=0,nalms-1
           l = nint(alms_in(i,1))
           m = nint(alms_in(i,2))
           if (m < 0) then
              write(*,'(a,i6,a,i6)') '     the a_lm with m<0 are ignored, l=',l,', m=',m
           else
  ! MR: the constraining a_lm in the input file could have entries with
  !     l>nlmax. If they were entered into alm_TGC, this would cause memory
  !     corruption, so this must be checked.
              if (l <= nlmax) then
                 alm_TGC(j,l,m) = &
!                &   CMPLX(alms_in(i,3)+randgauss_boxmuller(iseed)*alms_in(i,5), &
!                &         alms_in(i,4)+randgauss_boxmuller(iseed)*alms_in(i,6), kind=KALM) * beam_const(l,j)
               &   CMPLX(alms_in(i,3)+rand_gauss(rng_handle)*alms_in(i,5), &
               &         alms_in(i,4)+rand_gauss(rng_handle)*alms_in(i,6), kind=KALM) * beam_const(l,j)
              endif
           endif
        enddo
     enddo

     deallocate(alms_in)
     deallocate(beam_const)
  endif

  if (output_alms) then

     PRINT *,"      "//code//"> outputting alms "

     ! write each component (T,G,C) in a different extension of the same file
     do i = 1, n_pols
        header = ""
        ! put inherited information immediatly, so that keyword values can be updated later on
        ! by current code values
        call add_card(header,"COMMENT","****************************************************************")
        call merge_headers(header_PS, header) ! insert header_PS in header at this point
        call add_card(header,"COMMENT","****************************************************************")

        ! start putting information relative to this code and run
        call write_minimal_header(header, 'alm', append=.true., &
             creator = CODE, version = VERSION, &
             nlmax = nlmax, nmmax = nmmax, randseed = ioriginseed, &
             units = units_alm(i), polar = polarisation)
        if (i == 1) then
           call add_card(header,"EXTNAME","'SIMULATED a_lms (TEMPERATURE)'", update = .true.)
        elseif (i == 2) then
           call add_card(header,"EXTNAME","'SIMULATED a_lms (GRAD / ELECTRIC component)'", update = .true.)
        elseif (i == 3) then
           call add_card(header,"EXTNAME","'SIMULATED a_lms (CURL / MAGNETIC component)'", update = .true.)
        endif
        if (input_alms) then
           call add_card(header,"HISTORY","File constraining alms:")
           call add_card(header,"HISTORY",trim(infile_alms))
           if (apply_windows) then
              call add_card(header,"HISTORY","These alms are multiplied by pixel and beam window functions")
           else
              call add_card(header,"HISTORY","These alms are used as they are.")
           endif
        endif
        if (input_cl) then
           call add_card(header,"HISTORY","File with input C(l):")
           call add_card(header,"HISTORY",trim(infile))
           call add_card(header,"HISTORY","These alms are multiplied by pixel and beam window functions")
        endif
        if (input_cl .or. (input_alms.and.apply_windows)) then
           call add_card(header,"NSIDE"   ,nsmax,   "Resolution parameter for HEALPIX")
           if (trim(beam_file) == '') then
              call add_card(header,"FWHM"    ,fwhm_deg   ," [deg] FWHM of gaussian symmetric beam")
           else
              call add_card(header,"BEAM_LEG",trim(beam_file), &
                   & "File containing Legendre transform of symmetric beam")
           endif
        endif
        call add_card(header,"HISTORY")
        nlheader = SIZE(header)
        call dump_alms (outfile_alms,alm_TGC(i,0:nlmax,0:nmmax),nlmax,header,nlheader,i-1_i4b)
     enddo
  endif

  !   skip to the end if no map is asked for
  if (.not. do_map) goto 1000
  !-----------------------------------------------------------------------
  !                    precomputed plms
  !-----------------------------------------------------------------------

  if (n_plm.ne.0) then

     allocate(plm(0:n_plm-1,1:plm_pol),stat = status)
     call assert_alloc(status,code,"plm")

     if (plm_pol.eq.1) then
        PRINT*,"      "//code//"> Reading precomputed scalar P_lms "
     else
        PRINT*,"      "//code//"> Reading precomputed tensor P_lms "
     endif
!!     call read_dbintab(infile_plm,plm,n_plm,plm_pol,nullval,anynull=bad)
     call read_bintab(infile_plm, plm, n_plm, plm_pol, nullval, anynull=bad)
     if (bad) call fatal_error("Missing points in P_lm-file!")

  endif

  !-----------------------------------------------------------------------
  !                           a_lm to map
  !-----------------------------------------------------------------------

  PRINT *,"      "//code//"> Generating sky map(s) "

  zbounds=(/-1.0_DP, 1.0_DP/) ! full sky  = default
  select case (plm_pol + polar*10 + deriv*100)
  case(0) ! temperature only
     call alm2map(nsmax,nlmax,nmmax,alm_TGC,map_TQU(:,1), &
          &       zbounds=zbounds)
  case(1) ! temperature only with precomputed Ylm
     call alm2map(nsmax,nlmax,nmmax,alm_TGC,map_TQU(:,1),plm=plm(:,1))
  case(10) ! T+P
     call alm2map(nsmax,nlmax,nmmax,alm_TGC,map_TQU, &
                  zbounds=zbounds)
  case(11) ! T+P with precomputed Ylm_T
     call alm2map(nsmax,nlmax,nmmax,alm_TGC,map_TQU,plm=plm(:,1))
   case(13)  ! T+P with precomputed Ylm
     call alm2map(nsmax,nlmax,nmmax,alm_TGC,map_TQU,plm=plm)
  case(100)
     call alm2map_der(nsmax,nlmax,nmmax,alm_TGC,map_TQU(:,1),&
                      der1=map_TQU(:,2:3), zbounds=zbounds)
  case(200)
     call alm2map_der(nsmax,nlmax,nmmax,alm_TGC,map_TQU(:,1),&
          &           der1=map_TQU(:,2:3),der2=map_TQU(:,4:6), zbounds=zbounds)
  case(110)
     call alm2map_der(nsmax,nlmax,nmmax,alm_TGC,map_TQU(:,1:3),&
                      der1=map_TQU(:,4:9), zbounds=zbounds)
  case(210)
     call alm2map_der(nsmax,nlmax,nmmax,alm_TGC,map_TQU(:,1:3),&
          &           der1=map_TQU(:,4:9),der2=map_TQU(:,10:18), zbounds=zbounds)
  case default
     print*,'Not valid case'
     print*,'Polarisation:', polarisation
     print*,'Derivatives: ', deriv
     print*,'P_lm pol:    ', plm_pol
     call fatal_error(code)
  end select

  if (n_plm.ne.0) then
     deallocate(plm)
  endif

  !-----------------------------------------------------------------------
  !                      write the map to FITS file
  !-----------------------------------------------------------------------

  nlheader = SIZE(header)
  header = ""
  PRINT *,"      "//code//"> Writing sky map to FITS file "
  ! put inherited information immediatly, so that keyword values can be updated later on
  ! by current code values
  call add_card(header,"COMMENT","****************************************************************")
  call merge_headers(header_PS, header) ! insert header_PS in header at this point
  call add_card(header,"COMMENT","****************************************************************")
  ! start putting information relative to this code and run
  call write_minimal_header(header, 'map', append=.true., &
       nside = nsmax, ordering = 'RING', coordsys = coordsys, &
       fwhm_degree = fwhm_arcmin / 60.d0, &
       beam_leg = trim(beam_file), &
       polar = polarisation, &
       deriv = deriv, &
       creator = CODE, version = VERSION, &
       nlmax = nlmax, &
       randseed = ioriginseed, &
       units = units_map(1) )
  ! update 'EXTNAME' in place
  call add_card(header,"EXTNAME","'SIMULATED MAP'", update=.true.)
  call add_card(header,"COMMENT","*************************************")
  if (input_alms) then
     call add_card(header,"HISTORY","File constraining alms:")
     call add_card(header,"HISTORY",trim(infile_alms))
     if (apply_windows) then
        call add_card(header,"HISTORY","These alms are multiplied by pixel and beam window functions")
     else
        call add_card(header,"HISTORY","These alms are used as they are.")
     endif
  endif


  ! call output_map(map_TQU(0:npixtot-1,1:n_maps), header, outfile) !not on Sun
  call write_bintab(map_TQU, npixtot, n_maps, header, nlheader, outfile)
  
  !-----------------------------------------------------------------------
  !                      deallocate memory for arrays
  !-----------------------------------------------------------------------

  DEALLOCATE( map_TQU )
1000 continue ! no map done
  DEALLOCATE( alm_TGC )
  DEALLOCATE( units_alm, units_map )

  !-----------------------------------------------------------------------
  !                      output and report card
  !-----------------------------------------------------------------------
  call wall_clock_time(time1)
  call cpu_time(ptime1)
  clock_time = time1 - time0
  ptime      = ptime1 - ptime0

  write(*,9000) " "
  write(*,9000) " Report Card for "//code//" simulation run"
  write(*,9000) "----------------------------------------"
  write(*,9000) " "
  if (.not. polarisation) then
     chline = "      Temperature alone"
     if (deriv == 1) chline = "      Temperature + 1st derivatives"
     if (deriv == 2) chline = "     Temperature + 1st & 2nd derivatives"
  else
     chline = "    Temperature + Polarisation"
     if (deriv == 1) chline = "      T+P + 1st derivatives"
     if (deriv == 2) chline = "     T+P + 1st & 2nd derivatives"
  endif
  write(*,9000) chline
  write(*,9000) " "
  write(*,9000) " Input power spectrum : "//TRIM(infile)
  if (input_alms) then
     write(*,9000) " Constraining a_lms : "//TRIM(infile_alms)
  endif
  write(*,9005) " Multipole range      : 0 < l <= ", nlmax
  write(*,9010) " Number of pixels     : ", npixtot
  write(*,9020) " Pixel size in arcmin : ", pix_size_arcmin
  write(*,9010) " Initial random # seed: ", ioriginseed
  if (trim(beam_file) == '') then
     write(*,9020) " Gauss. FWHM in arcmin: ", fwhm_arcmin
  else
     write(*,9000) " Beam file: "//trim(beam_file)
  endif
  if (do_map) then
     write(*,9000) " Output map           : "//TRIM(outfile)
  else
     write(*,9000) " no map done!"
  endif
  if (output_alms) then
     write(*,9000) " Output a_lms in      : "//TRIM(outfile_alms)
  endif
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

