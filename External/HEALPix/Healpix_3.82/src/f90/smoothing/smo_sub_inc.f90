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
! Fast Smoothing of a full sky map in HEALPIX pixelisation.
! Written and developed by E. Hivon (hivon@iap.fr) and K. Gorski
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
  !     this code uses the FITSIO library that can be found at
  !     http://heasarc.gsfc.nasa.gov/docs/software/fitsio/fitsio.html
  !
  !  RELATED LITTERATURE:
  !     about HEALPIX   : see Gorski et al, 2005, ApJ, 622, 759
  !     about this code : see Hivon & Gorski, 1997, in preparation
  !
  !  HISTORY:
  !     April-October 1997, Eric Hivon, TAC (f77)
  !               Jan 1998 : translation in Fortran 90
  !               Dec 2001 : version 1.2, EH, IPAC
  !                 implementation of non-interactive input
  !                 extension to polarised maps
  !               Sep 2002 : implement new parser
  !               Feb 2003 : output correct ORDERING keyword
  !               May 2007 : calls map2alm_iterative
  !               Jun 2010 : supports large maps
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
  integer(I4B),             parameter :: KALMC = KIND((1.0_KMAP, 1.0_KMAP)) ! precision of alm arrays
  !-------------------------------------------------------------------------

  integer(I4B) :: nsmax, nlmax, nmmax

  complex(kind=KALMC), DIMENSION(:,:,:), ALLOCATABLE :: alm_TGC
  complex(kind=KALMC), DIMENSION(:,:,:), ALLOCATABLE :: alm_TGC2
  real(kind=KMAP),     DIMENSION(:,:),   ALLOCATABLE :: map_TQU
  real(kind=KMAP),     DIMENSION(:,:),   ALLOCATABLE :: map_TQU2
  real(kind=DP),       DIMENSION(:,:),   ALLOCATABLE :: w8ring_TQU
  real(kind=DP),       DIMENSION(:,:),   ALLOCATABLE :: dw8
  real(kind=DP),       DIMENSION(:,:),   ALLOCATABLE :: plm
!  real(kind=KMAP),     dimension(:,:),   ALLOCATABLE :: pw8list
  real(kind=KMAP),     dimension(:,:), target, ALLOCATABLE :: pw8map
  real(kind=KMAP),     dimension(:,:), pointer       :: pmask
  integer(kind=I4B) :: status

  integer(kind=I8B) :: npixtot, n_plm, ipix, npixtot_tmp !, n_pw8_read
  integer(kind=I4B) :: mlpol
  integer(kind=I4B) :: plm_nside, plm_lmax, plm_pol, plm_mmax
  !integer(kind=I4B) :: nmaps, type, polar_fits
  integer(kind=I4B) :: i, j
  integer(kind=I4B) :: nlheader
  integer(kind=I4B) :: ordering, order_type
  integer(kind=I4B) :: iter, iter_order, polar, n_pols
  integer(kind=I4B) :: simul_type
  integer(kind=I4B) :: n_rings, n_rings_read, won, nmw8

  real(kind=KMAP)               :: fwhm_arcmin, fwhm_deg
  real(kind=KMAP)               :: fmissval
  real(kind=SP),      parameter :: fbad_value = -1.6375e30_sp
  real(kind=SP)                 :: pix_size_arcmin
  real(kind=DP)                 :: theta_cut_deg, cos_theta_cut
  real(kind=DP)                 :: nullval
  real(kind=DP), dimension(1:2) :: zbounds
  real(kind=DP)                 :: fsky
  real(kind=SP)                 :: clock_time, time0, time1
  real(kind=SP)                 :: ptime, ptime0, ptime1

  character(LEN=filenamelen)          :: infile
  character(LEN=filenamelen)          :: outfile
!  character(LEN=filenamelen)          :: parafile = ''
  character(LEN=filenamelen)          :: infile_plm
  character(LEN=filenamelen)          :: infile_w8
  character(LEN=filenamelen)          :: w8name
  character(LEN=filenamelen)          :: def_dir, def_file
  character(LEN=filenamelen)          :: usr_dir, usr_file
  character(LEN=filenamelen)          :: final_file
  character(LEN=filenamelen)          :: healpixtestdir
  character(LEN=filenamelen)          :: beam_file
  character(len=filenamelen)          :: description
  character(len=100)                  :: chline, chline1, strw8file
  character(LEN=80), DIMENSION(1:180) :: header
  LOGICAL(kind=LGT) :: bad, ok, polarisation
!   character(LEN=*), PARAMETER :: code = 'SMOOTHING'
  character(len=*), parameter :: VERSION = HEALPIX_VERSION
  integer(kind=i4b), parameter :: nm_max = 5
  character(len=80), dimension(1:nm_max)   :: units_map
  character(len=20)                        :: coordsys

  type(paramfile_handle) :: handle

  real(kind=DP), dimension(0:3) :: mono_dip
  real(kind=DP)                 :: theta_dip, phi_dip, long_dip, lat_dip
  integer(kind=I4B)             :: lowlreg

  !-----------------------------------------------------------------------
  !                    get input parameters and create arrays
  !-----------------------------------------------------------------------

  call wall_clock_time(time0)
  call cpu_time(ptime0)
  !     --- read parameters interactively if no command-line arguments
  !     --- are given, otherwise interpret first command-line argument
  !     --- as parameter file name and read parameters from it:

  !     --- announces program, default banner ***
  print *, " "
  print *,'                        '//code//' '//version
  write(*,'(a)') '         *** Smoothing of a Temperature (+Polarisation) map ***    '
  if (KMAP == SP) print*,'                Single precision outputs'
  if (KMAP == DP) print*,'                Double precision outputs'
  print *, ''
  handle = parse_init(parafile)

  !     --- choose temp. only or  temp. + pol. ---
  description = concatnl( &
       & " Do you want to analyse", &
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

  !     --- check that this map is usable for this analysis
  call check_input_map(code, infile, polarisation) ! this may change polarisation from T to F

  !     --- finds out the pixel number of the map and its ordering ---
  ! npixtot = getsize_fits(infile, nmaps = nmaps, ordering=ordering, nside=nsmax, &
  ! &    mlpol=mlpol, type = type, polarisation = polar_fits, coordsys = coordsys)

  !    --- get other relevant information ---
  npixtot_tmp = getsize_fits(infile, ordering=ordering, nside=nsmax, &
       mlpol=mlpol, coordsys = coordsys)

!   if (nsmax<=0) then
!      print*,'Keyword NSIDE not found in FITS header!'
!      call fatal_error(code)
!   endif
!   if (type == 3) npixtot = nside2npix(nsmax) ! cut sky input data set
!   if (nsmax/=npix2nside(npixtot)) then
!      print*,"FITS header keyword NSIDE does not correspond"
!      print*,"to the size of the map!"
!      call fatal_error(code)
!   endif

!   if (polarisation .and. (nmaps >=3) .and. polar_fits == -1) then
!      print*,"The input fits file MAY NOT contain polarisation data."
!      print*,"Proceed at your own risk"
!   endif

!   if (polarisation .and. (nmaps<3 .or. polar_fits ==0)) then
!      print *,"The file does NOT contain polarisation maps"
!      print *,"only the temperature field will be smoothed"
!      polarisation = .false.
!   endif

!   !     --- check ordering scheme ---
!   if ((ordering/=1).and.(ordering/=2)) then
!      print*,'The ordering scheme of the map must be RING or NESTED.'
!      print*,'No ordering specification is given in the FITS-header!'
!      call fatal_error(code)
!   endif

  order_type = ordering

  !     --- gets the L range for the analysis ---
  WRITE(chline1,"(a,i5)") " The map has Nside = ",nsmax
  if (mlpol<=0) WRITE(chline,"(a,i5,a)") "We recommend: (0 <= l <= l_max <= ",3*nsmax-1,")"
  if (mlpol> 0) WRITE(chline,"(a,i5,a)") "We recommend: (0 <= l <= l_max <= ",mlpol,")"
  description = concatnl(&
       & "", &
       & chline1, &
       & " Enter the maximum l range (l_max) for the simulation. ", &
       & chline )
  nlmax = parse_int(handle, 'nlmax', default=2*nsmax, descr=description)

  ! ----------  default parameters ----------------
  nmmax   = nlmax
  npixtot = nside2npix(nsmax)
  pix_size_arcmin = 360.*60./SQRT(npixtot*PI)

  !     --- ask about monopole/dipole removal ---
  description = concatnl(&
       & "", &
       & " To improve the temperature map analysis (specially on a cut sky) ", &
       & " you have the option of first regressing out the best fit monopole and dipole", &
       & "   NOTE : the regression is made on valid (unflagged) pixels, ", &
       & " out of the symmetric cut (if any)", &
       & "  Do you want to : ", &
       & " 0) Do the analysis on the raw map", &
       & " 1) Remove the monopole (l=0) first ", &
       & " 2) Remove the monopole and dipole (l=1) first")
  lowlreg = parse_int(handle, 'regression', vmin=0, vmax=2, default=0, descr=description)


  !     --- gets the cut applied to the map ---
  description = concatnl(&
       & "", &
       & " Enter the symmetric cut around the equator in DEGREES : ", &
       & " (One ignores data within |b| < b_cut)     0 <= b_cut = ")
  theta_cut_deg = parse_double(handle, 'theta_cut_deg', &
       &                       vmin=0.0_dp, default=0.0_dp, descr=description)
  cos_theta_cut = SIN(theta_cut_deg/180.d0*PI) !counted from equator instead of from pole

  zbounds = (/ cos_theta_cut , -cos_theta_cut /)
  if (theta_cut_deg<1e-4) zbounds = (/ -1.0_dp, 1.0_dp /) !keep all sphere if cut not set
  fsky = (zbounds(2)-zbounds(1))/2.0_dp
  if (fsky <= 0.0_dp) fsky = 1.0_dp + fsky
  write(*,"(a,f6.1,a)") "One keeps ",100.*fsky," % of the original map"

  !     --- gets the order of iteration ---
  description = concatnl(&
       & "", &
       & " Do you want : ", &
       & " 0) a standard analysis", &
       & " 1,2,3,4....) an iterative analysis", &
       & " (enter order of iteration, 3rd order is usually optimal)")
  iter_order=parse_int(handle, 'iter_order', vmin=0, default=0, descr=description)
  if (.not. handle%interactive) then
     select case (iter_order)
     case (0)
        print*," Standard analysis"
     case default
        print*," Iterative analysis"
     end select
  end if

  ! ------------------- precomputed plms  ---------------------
40 continue
  n_plm   = 0
  plm_pol = 0
  description = concatnl(&
       & "", &
       & " Enter name of file with precomputed p_lms :", &
       & "(eg, plm.fits) (if '' is typed, no file is read)" )
  infile_plm = parse_string(handle, 'plmfile', default="''", descr=description, filestatus='old')

  if (trim(infile_plm) /= "") then
     ! check that the plm file matches the current simulation (for nside,lmax,mmax)
     call read_par(infile_plm, plm_nside, plm_lmax, plm_pol, mmax=plm_mmax)

     if ((plm_pol /= 1).and.(.not.polarisation)) plm_pol=1

     if (plm_nside /= nsmax .or. plm_lmax /= nlmax .or. plm_mmax /= nmmax) then
        print *," "//code//"> Plm file does not match map to smooth in Nside, Lmax or Mmax ! "
        print*,'nside (plm file, map) = ',plm_nside,nsmax
        print*,'lmax                  = ',plm_lmax, nlmax
        print*,'mmax                  = ',plm_mmax, nmmax
        if (handle%interactive) goto 40
        call fatal_error(code)
     endif
     n_plm = (nmmax + 1_i8b)*(2*nlmax + 2_i8b -nmmax)*nsmax
     print *," "
  endif

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
  if (beam_file /= '') then
     fwhm_arcmin = 0.
     print*,'fwhm_arcmin is now : 0.'
  endif

  !     --- gets the output sky map filename ---
  description = concatnl(&
       & "", &
       & " Enter Output map file name (eg, test_smooth.fits) :", &
       & "  (or !test_smooth.fits to overwrite an existing file)" )
  outfile = parse_string(handle, "outfile", &
       default="!test_smooth.fits", descr=description, filestatus="new")

  !-----------------------------------------------------------------------
  !                      ask for weights
  !-----------------------------------------------------------------------
  description = concatnl(&
       & "", &
       & " Do you want to use quadrature weights in the analysis: ", &
       & " 0: no weight            ", &
       & " 1: ring-based weights   ", &
       & " 2: pixel-based weights ?")
  won = parse_int(handle, 'won', vmin=0, vmax=2, default=0, descr=description)

  infile_w8=""
  strw8file = ""

  if (won == 1 .or. won == 2) then

     ! default weight file name
     w8name = get_healpix_weight_file(nsmax, won)

     def_file = trim(w8name)
     def_dir  = get_healpix_data_dir()

     if (won == 1) strw8file = "ring-based  quadrature weight file"
     if (won == 2) strw8file = "pixel-based quadrature weight file"

22   continue
     final_file = ''
     ok = .false.
     ! if interactive, try default name in default directories first
     if (handle%interactive) ok = scan_directories(def_dir, def_file, final_file)
     if (.not. ok) then
        ! not found, ask the user
        description = concatnl("",&
             &        " Could not find "//trim(strw8file), &
             &        " Enter the directory where this file can be found:")
        usr_dir = parse_string(handle,'w8filedir',default='',descr=description)
        if (trim(usr_dir) == '') usr_dir = trim(def_dir)
        description = concatnl("",&
             &        " Enter the name of the "//trim(strw8file)//":")
        usr_file = parse_string(handle,'w8file',default=def_file,descr=description)
        ! look for new name in user provided or default directories
        ok   = scan_directories(usr_dir, usr_file, final_file)
        ! if still fails, crash or ask again if interactive
        if (.not. ok) then
           print*,' File not found:'//trim(usr_file)
           if (handle%interactive) goto 22
           call fatal_error(code)
        endif
     endif
     infile_w8 = final_file
  endif

  print *," "
  call parse_check_unused(handle, code=lcode)
  call parse_summarize(handle,code=lcode,prec=KMAP)
  call parse_finish(handle)
!  call parse_summarize(handle)
  call brag_openmp()

  !-----------------------------------------------------------------------
  !              allocate space for arrays
  !-----------------------------------------------------------------------
  polar = 0
  if (polarisation) polar = 1
  n_pols = 1 + 2*polar ! either 1 or 3

  ALLOCATE(map_TQU(0:12*nsmax**2-1,1:n_pols),stat = status)
  call assert_alloc(status,code,"map_TQU")

  ALLOCATE(w8ring_TQU(1:2*nsmax,1:n_pols),stat = status)
  call assert_alloc(status,code,"w8ring_TQU")

  ALLOCATE(dw8(1:2*nsmax,1:n_pols),stat = status)
  call assert_alloc(status,code,"dw8")

  if (n_plm/=0) then
     ALLOCATE(plm(0:n_plm-1,1:plm_pol),stat = status)
     call assert_alloc(status,code,"plm")
  endif
  !-----------------------------------------------------------------------
  !                      reads the map
  !-----------------------------------------------------------------------
  print *,"      "//code//"> Inputting original map "

  fmissval = 0.0
  if (lowlreg > 0) fmissval = fbad_value
  call input_map(infile, map_TQU(0:,1:n_pols), &
       &  npixtot, n_pols, fmissval=fmissval, units=units_map(1:n_pols))
  do i=1,n_pols
     units_map(i) = adjustl(units_map(i))
     if (trim(units_map(i)) == '') units_map(i) = 'unknown'
  enddo
 
  !-----------------------------------------------------------------------
  !                      remove dipole
  !-----------------------------------------------------------------------

  if (lowlreg > 0) then
     PRINT *,"      "//code//"> Remove monopole (and dipole) from Temperature map"
     call remove_dipole(nsmax, map_TQU(0:npixtot-1,1), ordering, lowlreg, &
          & mono_dip(0:), zbounds, fmissval=fmissval)

     if (fmissval /= 0.0) then
        do j=1,n_pols
           do ipix=0, npixtot-1
!               if (abs(map_TQU(i,j)/fmissval-1.0) < 1.e-6*fmissval) then
              if (abs(map_TQU(ipix,j)-fmissval) <= abs(1.e-6*fmissval)) then
                 map_TQU(ipix,j) = 0.0_sp
              endif
           enddo
        enddo
     endif

     write(unit=*,fmt="(a,g13.3,a)")  " Monopole = ",&
          & mono_dip(0)," "//trim(units_map(1))
     if (lowlreg > 1) then
        call vec2ang( mono_dip(1:3), theta_dip, phi_dip)
        write(unit=*,fmt="(a,g13.3,a)")  " Dipole   = ",&
             & sqrt(sum(mono_dip(1:3)**2))," "//trim(units_map(1))
        long_dip =      phi_dip  /PI*180.0
        lat_dip  = 90.0-theta_dip/PI*180.0
        write(unit=*,fmt="(a,g10.3,', ',g10.3,a)")  &
             & "(long.,lat.) = (",long_dip,lat_dip,") Deg"
     endif

  endif
  !-----------------------------------------------------------------------
  !                      input weights
  !-----------------------------------------------------------------------

  dw8=0.0_dp ! read as DP, even in SP in FITS file
  if (trim(infile_w8)/="") then
     if (won == 1) then
        ! ring based weights: treated as such in map2alm
        n_rings = 2 * nsmax
        n_rings_read = getsize_fits(infile_w8, nmaps=nmw8)

        if (n_rings_read /= n_rings) then
           print *," "
           print *,"wrong ring weight file:"//trim(infile_w8)
           call fatal_error(code)
        endif
        nmw8 = min(nmw8, n_pols)

        PRINT *,"      "//code//"> Inputting "//trim(strw8file)
        call input_map(infile_w8, dw8, n_rings, nmw8, fmissval=0.0_dp)
     endif

     if (won == 2) then 
        ! pixel based weights: turned into full sky map and treat as mask
        allocate(pw8map (0:npixtot-1,   1:1), stat=status)
        call assert_alloc(status,code,"pw8map")

        PRINT *,"      "//code//"> Inputting and unfolding "//trim(strw8file)

!         n_pw8_read = getsize_fits(infile_w8)
!         allocate(pw8list(0:n_pw8_read-1,1:1), stat=status)
!         call assert_alloc(status,code,"pw8list")
!         call input_map(infile_w8, pw8list, n_pw8_read, 1)
!         call unfold_weightslist(nsmax, won, pw8list(:,1), pw8map(:,1))
!         deallocate(pw8list)

        call unfold_weightsfile(infile_w8, pw8map(:,1))

        pmask => pw8map(:,1:1)
        !print*,pmask(1:12,1)
        !print*,pmask(npixtot-11:npixtot,1)
     endif
  endif
  w8ring_TQU = 1.d0 + dw8
  if (won /= 2) allocate(pmask(0:1,1))

  !-----------------------------------------------------------------------
  !                 reorder maps to RING if necessary
  !-----------------------------------------------------------------------
  if (order_type == 2) then
     print *,"      "//code//"> Convert Nest -> Ring "
     call convert_nest2ring (nsmax, map_TQU)
  endif

  !-----------------------------------------------------------------------
  !                    precomputed plms
  !-----------------------------------------------------------------------

  if (n_plm/=0) then
     if (plm_pol.eq.1) then
        print*,"      "//code//"> Reading precomputed scalar P_lms "
     else
        print*,"      "//code//"> Reading precomputed tensor P_lms "
     endif
     !call read_dbintab(infile_plm,plm,n_plm,plm_pol,nullval,anynull=bad) ! replaced 2010-11-24
     call read_bintab(infile_plm, plm, n_plm, plm_pol, nullval, anynull=bad)
     if (bad) call fatal_error ("Missing points in P_lm-file!")
  endif

  !-----------------------------------------------------------------------
  !                        map to alm
  !-----------------------------------------------------------------------

  ALLOCATE(alm_TGC(1:n_pols, 0:nlmax, 0:nmmax),stat = status)
  call assert_alloc(status,code,"alm_TGC")

  PRINT *,"      "//code//"> Analyse map "

  if (n_plm/=0) then
     call map2alm_iterative(nsmax, nlmax, nmmax, iter_order, map_TQU, alm_TGC, &
          &                 zbounds=zbounds, w8ring=w8ring_TQU, plm=plm, mask=pmask)
  else
     call map2alm_iterative(nsmax, nlmax, nmmax, iter_order, map_TQU, alm_TGC, &
          &                 zbounds=zbounds, w8ring=w8ring_TQU, mask=pmask)
  endif

  !-----------------------------------------------------------------------
  !                    smoothing of the alm
  !-----------------------------------------------------------------------
  print *,'      '//code//'> Alter alm '
  call alter_alm(nsmax, nlmax, nmmax, fwhm_arcmin, alm_TGC(1:1+2*polar,0:nlmax,0:nmmax), beam_file)

  !-----------------------------------------------------------------------
  !                           a_lm to map
  !-----------------------------------------------------------------------

  print *,'      '//code//'> Generating sky map(s) '

  select case (plm_pol+polar*10)
  case(0)
     call alm2map(nsmax,nlmax,nmmax,alm_TGC,map_TQU(:,1))
  case(1)
     call alm2map(nsmax,nlmax,nmmax,alm_TGC,map_TQU(:,1),plm=plm(:,1))
  case(10)
     call alm2map(nsmax,nlmax,nmmax,alm_TGC,map_TQU)
  case(11)
     call alm2map(nsmax,nlmax,nmmax,alm_TGC,map_TQU,plm=plm(:,1))
  case(13)
     call alm2map(nsmax,nlmax,nmmax,alm_TGC,map_TQU,plm=plm)
  end select

  !-----------------------------------------------------------------------
  !                 convert back to NEST if necessary
  !-----------------------------------------------------------------------
  if (order_type == 2) then
     print *,"      "//code//"> Convert Ring -> Nest "
     call convert_ring2nest (nsmax, map_TQU)
  endif

  !-----------------------------------------------------------------------
  !                        generates header
  !-----------------------------------------------------------------------
  nlheader = SIZE(header)
  fwhm_deg = fwhm_arcmin/60.
  print *,'      '//code//'> Writing smoothed map to FITS file '

  call write_minimal_header(header, 'map', &
       order = order_type, nside = nsmax, coordsys = coordsys, &
       creator = CODE, version = VERSION, nlmax = nlmax, &
       beam_leg = trim(beam_file), fwhm_degree = fwhm_deg * 1.d0, &
       polar = polarisation, units = units_map(1) )

  call add_card(header,'EXTNAME','SMOOTHED DATA', update = .true.)
  call add_card(header)
  call add_card(header,'COMMENT','************************* Input data ************************* ')
  call add_card(header,'COMMENT','Input Map in '//TRIM(infile))
  call add_card(header,'COMMENT','*************************************************************** ')
  call add_card(header,"NITERALM",iter_order,   " Number of iterations for a_lms extraction")
  call add_card(header,"LREGRESS",lowlreg,      " Regression of low multipoles (0/1/2)")
  call add_card(header,"QUADSCHM",won,          " Quadrature scheme (0/1/2)")
  call add_card(header) ! blank line
  !-----------------------------------------------------------------------
  !                      write the map to FITS file
  !-----------------------------------------------------------------------
  !call output_map(map_TQU(0:npixtot-1,1:n_pols), header, outfile) !not on Sun
  call write_bintab(map_TQU, npixtot, n_pols, header, nlheader, outfile)

  !-----------------------------------------------------------------------
  !                      deallocate memory for arrays
  !-----------------------------------------------------------------------
  DEALLOCATE(map_TQU)
  DEALLOCATE(w8ring_TQU)
  DEALLOCATE(dw8)

  !-----------------------------------------------------------------------
  !                      output and report card
  !-----------------------------------------------------------------------
  call wall_clock_time(time1)
  call cpu_time(ptime1)
  clock_time = time1 - time0
  ptime      = ptime1 - ptime0

  write(*,9000) " "
  write(*,9000) " Report Card for "//code//" run"
  write(*,9000) " -----------------------------"
  write(*,9000) " "
  if (.not. polarisation) then
     write(*,9000) "      Temperature alone"
  else
     write(*,9000) "    Temperature + Polarisation"
  endif
  write(*,9000) " "

  write(*,9000) " Input map              : "//TRIM(infile)
  write(*,9006) " Multipole range        : ",lowlreg ," <= l <= ", nlmax
  write(*,9010) " Number of pixels       : ", npixtot
  write(*,9020) " Pixel size in arcmin   : ", pix_size_arcmin
  if (trim(beam_file) == '') then
     write(*,9020) " Gauss. FWHM in arcmin  : ", fwhm_arcmin
  else
     write(*,9000) " Beam file  : "//trim(beam_file)
  endif
  write(*,9000) " Output map             : "//TRIM(outfile)
  if (won >0) write(*,9000) " Quadrature weight file : "//trim(infile_w8)
  write(*,9030) " Clock and CPU time [s] : ", clock_time, ptime

  !-----------------------------------------------------------------------
  !                       end of routine
  !-----------------------------------------------------------------------
  write(*,9000) " "
  write(*,9000) " "//code//"> normal completion"

9000 format(a)
9005 format(a,i8)
9006 format(a,i2,a,i8)
9010 format(a,i16)
9020 format(a,g20.5)
9030 format(a,f11.2,f11.2)
