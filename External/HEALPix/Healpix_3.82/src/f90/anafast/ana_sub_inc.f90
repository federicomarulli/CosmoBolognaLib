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
! Computes the spherical harmonics coefficients (a_lm)
! and the power spectrum from a full sky map in HEALPIX pixelisation
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
  !     this code uses the CFITSIO library that can be found at
  !     http://heasarc.gsfc.nasa.gov/docs/software/fitsio/
  !
  !  RELATED LITTERATURE:
  !     about HEALPIX   : see Gorski et al, 2005, ApJ, 622, 759
  !     about this code : see Hivon & Gorski, in preparation
  !
  !  HISTORY:
  !     April-October 1997, Eric Hivon, TAC (f77)
  !               Jan 1998 : translation in Fortran 90
  !               Aug 1998 : version 1.0.0
  !               Apr 1999 : precomp. plms, higher order iterations and alm dumping, FKH
  !               Aug 2000 : CMBFAST convention for polarisation
  !               Sep 2002 : implement new parser
  !               May 2007 : can cross-correlate 2 maps, introduced map2alm_iterative
  !               Jun 2010 : supports large maps
  !               May 2012:  pro/down-grading of mask now works
  !
  !  FEEDBACK:
  !     for any questions : hivon@iap.fr
  !
  !=======================================================================
  !     version 2.3
  !=======================================================================
  ! this file can not be compiled on its own.
  ! It must be inserted into the file anafast.f90 by the command  include
  !
  integer(I4B),             parameter :: KALMC = KIND((1.0_KMAP, 1.0_KMAP)) ! precision of alm arrays
  integer(I4B),             parameter :: KCL   = KMAP  ! precision of c(l) arrays
  !-------------------------------------------------------------------------

  integer(I4B) :: nsmax, nlmax, nmmax, nsmax2

  complex(kind=KALMC), DIMENSION(:,:,:), ALLOCATABLE :: alm_TGC, alm_TGC2
  real(kind=KMAP),     DIMENSION(:,:),   ALLOCATABLE :: map_TQU, map_TQU2
  real(kind=KMAP),     DIMENSION(:,:), target, ALLOCATABLE :: mask, mask_tmp
  real(kind=KMAP),     dimension(:,:), pointer :: pmask, pmask_tmp
  real(kind=KCL),      DIMENSION(:,:),   ALLOCATABLE :: clout
  real(kind=KCL),      DIMENSION(:),     ALLOCATABLE :: fact_norm
  real(kind=DP),       DIMENSION(:,:),   ALLOCATABLE :: w8ring_TQU
  real(kind=DP),       DIMENSION(:,:),   ALLOCATABLE :: dw8
  real(kind=DP),       DIMENSION(:,:),   ALLOCATABLE :: plm
!  real(kind=KMAP),     dimension(:,:),   ALLOCATABLE :: pw8list
  real(kind=KMAP),     dimension(:,:), target, ALLOCATABLE :: pw8map
  integer(kind=I4B) :: status

  integer(kind=I4B) :: mlpol, mlpol2
  integer(kind=I4B), dimension(1:2) :: order_map
  integer(kind=I4B) :: l, plm_nside, plm_lmax, plm_pol, plm_mmax
  integer(kind=I8B) :: npixtot, npixtot2, ipix, n_plm !, n_pw8_read
  integer(kind=I4B) :: i, j
  integer(kind=I4B) :: iter_order, polar, n_pols, n_maps
  integer(kind=I4B) :: simul_type
  integer(kind=I4B) :: n_rings, n_rings_read, won, nmw8

  real(kind=KMAP)               :: fmissval, fmiss_mask
  real(kind=DP)                 :: pix_size_arcmin
  real(kind=DP)                 :: theta_cut_deg, cos_theta_cut
  real(kind=DP)                 :: nullval
  real(kind=DP), dimension(1:2) :: zbounds
  real(kind=DP)                 :: fsky
  real(kind=SP)                 :: clock_time, time0, time1
  real(kind=SP)                 :: ptime, ptime0, ptime1
  real(kind=DP)                 :: min_mask, max_mask

  character(len=80), DIMENSION(1:120) :: header
  integer(kind=I4B) :: ncl, nlheader

  character(len=FILENAMELEN), dimension(2) :: infile
  character(len=FILENAMELEN)          :: outfile
!  character(len=FILENAMELEN)          :: parafile = ''
  character(len=FILENAMELEN), dimension(2) :: outfile_alms
  character(len=FILENAMELEN)          :: infile_plm
  character(len=FILENAMELEN)          :: infile_w8
  character(len=FILENAMELEN)          :: w8name
  character(len=FILENAMELEN)          :: def_dir, def_file
  character(len=FILENAMELEN)          :: usr_dir, usr_file
  character(len=FILENAMELEN)          :: final_file
  character(len=FILENAMELEN)          :: healpixtestdir
  character(len=FILENAMELEN)          :: description
  character(len=FILENAMELEN)          :: maskfile
  character(len=100)                  :: chline, chline1, strw8file
  character(len=20)                   :: coordsys, coordsys2
!  character(len=5)                    :: sstr
  LOGICAL(kind=LGT) :: bad, ok, polarisation, bcoupling, do_mask, twomaps, asym_cl
!   character(len=*), PARAMETER :: code = "ANAFAST"
  character(len=*), parameter :: VERSION = HEALPIX_VERSION
  character(len=80), dimension(1:3,1:2)   :: units_map
  character(len=80), dimension(1:3)       :: units_pow

  type(paramfile_handle) :: handle

  real(kind=DP), dimension(0:3,1:2) :: mono_dip
  real(kind=DP)                 :: theta_dip, phi_dip, long_dip, lat_dip, ampl_dip
  integer(kind=I4B)             :: lowlreg

  integer(kind=I4B)        :: npix_mask, order_mask, nside_mask, nmasks, i_mask
  logical(LGT)             :: pessim_mask
  !=========================================================================

  call wall_clock_time(time0)
  call cpu_time(ptime0)
  !-----------------------------------------------------------------------
  !                    get input parameters
  !-----------------------------------------------------------------------
  !     --- read parameters interactively if no command-line arguments
  !     --- are given, otherwise interpret first command-line argument
  !     --- as parameter file name and read parameters from it:

  !     --- announces program, default banner ***
  PRINT *, " "
  PRINT *, "                "//code//" "//version
  write(*,'(a)') " *** Analysis and power spectrum calculation for a Temperature/Polarisation map ***"
  if (KMAP == SP) print*,'                Single precision outputs'
  if (KMAP == DP) print*,'                Double precision outputs'
  PRINT *, " "
  handle = parse_init(parafile)

  !     --- choose temp. only or  temp. + pol. ---
  description = concatnl( &
       & " Do you want to analyse", &
       & " 1) a Temperature only map", &
       & " 2) Temperature + polarisation maps")
  simul_type = parse_int(handle, 'simul_type', default=1, vmin=1, vmax=2, descr=description)
  if (simul_type .eq. 1) polarisation = .false.
  if (simul_type .eq. 2) polarisation = .true.
  bcoupling = polarisation ! <<<<
!!  print*, " Bcoupling :  ",bcoupling
  order_map = -1

  !     --- gets the file name for the map(s) ---
  chline = "''"
  healpixtestdir = get_healpix_test_dir()
  if (trim(healpixtestdir)/='') chline = trim(healpixtestdir)//'/map.fits'
  description = concatnl( &
       & "", &
       & " Enter input file name 1 (Map FITS file): ")
  infile(1) = parse_string(handle, 'infile',default=chline, descr=description, filestatus='old')
  n_maps = 1
  description = concatnl( &
       & "", &
       & " Enter input file name 2 (Map FITS file): ")
  infile(2) = parse_string(handle, 'infile2',default="''", descr=description, filestatus='old')
  twomaps = (trim(infile(2)) /= "")
  if (twomaps) n_maps = 2

  !     --- check that this map is usable for this analysis
  call check_input_map(code, infile(1), polarisation) ! this may change polarisation from T to F

  !     --- finds out the pixel number of the map and its ordering and coordinates ---
  npixtot = getsize_fits(infile(1), ordering=order_map(1), &
       &                 nside=nsmax, mlpol=mlpol, coordsys=coordsys)
  if (trim(coordsys) == '') coordsys = 'unknown'

  if (twomaps) then
     call check_input_map(code, infile(2), polarisation)
     npixtot2 = getsize_fits(infile(2), ordering=order_map(2), &
          &                  nside=nsmax2, mlpol=mlpol2, coordsys=coordsys2)
     if (trim(coordsys2) == '') coordsys2 = 'unknown'
     mlpol = max(mlpol, mlpol2)
     if (nsmax2 /= nsmax .or. trim(coordsys2) /= trim(coordsys)) then
        print *," "//code//"> The 2 maps to be analyzed"
        print 9000, trim(infile(1)), trim(infile(2))
        print*,'do not match in resolution or coordinates'
        print*,nsmax,nsmax2
        print*,trim(coordsys)//' '//trim(coordsys2)
        call fatal_error(code)
     endif
  endif

  !     --- gets the L range for the analysis ---
  WRITE(chline1,"(a,i5)") " The map has Nside = ",nsmax
  if (mlpol<=0) WRITE(chline,"(a,i5,a)") "We recommend: (0 <= l <= l_max <= ",3*nsmax-1,")"
  if (mlpol> 0) WRITE(chline,"(a,i5,a)") "We recommend: (0 <= l <= l_max <= ",mlpol,")"
  description = concatnl(&
       & "", &
       & chline1, &
       & " Enter the maximum l range (l_max) for the analysis. ", &
       & chline )
  nlmax = parse_int(handle, 'nlmax', default=2*nsmax, descr=description)

  ! ----------  default parameters ----------------
  nmmax   = nlmax
  npixtot = nside2npix(nsmax)
  pix_size_arcmin = 360.*60./SQRT(npixtot*PI)
  !   -------------------------------------

  !     --- ask for pixel mask ---
  description = concatnl(&
       & "", &
       & " Enter FITS file containing pixel mask and/or weight you want to use: ", &
       & " the map(s) will be multiplied by this mask(s) before analysis ", &
       & " (If '', all valid pixels are used with their current value) ")
  maskfile = parse_string(handle, 'maskfile', default="''", descr=description, filestatus='old')
  do_mask = (trim(maskfile) /= "")

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
  write(*,"(a,f6.1,a)") "One keeps ",100.*fsky," % of the original map(s)"

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

  !          --- precomputed plms  ---
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
        print *," "//code//"> Plm file does not match map to analyze in Nside, Lmax or Mmax ! "
        print*,'nside (plm file, map) = ',plm_nside,nsmax
        print*,'lmax                  = ',plm_lmax, nlmax
        print*,'mmax                  = ',plm_mmax, nmmax
        if (handle%interactive) goto 40
        call fatal_error(code)
     endif
     n_plm = (nmmax + 1_i8b)*(2*nlmax + 2_i8b - nmmax)*nsmax
     print *," "
  endif

50 continue
  !     --- gets the output power spectrum filename ---
  if (twomaps) then
     chline =  " Enter Output CROSS Power Spectrum file name (eg, cl_out.fits) :"
  else
     chline =  " Enter Output Auto Power Spectrum file name (eg, cl_out.fits) :"
  endif
  description = concatnl(&
       & "", &
       & trim(chline), &
       & "  (or !cl_out.fits to overwrite an existing file)" )
  outfile = parse_string(handle, "outfile", &
       default="!cl_out.fits", descr=description, filestatus="new")

  !     --- gets the output alm-filename ---
  description = concatnl(&
       & "", &
       & " Enter output file name for alms(eg: alms.fits): ", &
       & " (or !alms.fits to overwrite an existing file): ", &
       & " (If '', the alms are not written to a file) ")
  outfile_alms(1) = parse_string(handle, 'outfile_alms', default="''", descr=description, filestatus='new')

  outfile_alms(2) = ''
  if (twomaps) then
     description = concatnl(&
          & "", &
          & " Enter output file name for alms of the second map")
     outfile_alms(2) = parse_string(handle, 'outfile_alms2', default="''", descr=description, filestatus='new')
  endif

  if (trim(outfile)=="" .and. trim(outfile_alms(1))=="" .and. trim(outfile_alms(2))=="") then
     print*,'ERROR: No output file provided'
     if (handle%interactive) goto 50
     call fatal_error(code)
  endif
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
        PRINT*," Standard analysis"
     case default
        PRINT*," Iterative analysis"
     end select
  end if

  PRINT *," "
  call parse_check_unused(handle, code=lcode)
  call parse_summarize(handle,code=lcode,prec=KMAP)
  call parse_finish(handle)
  call brag_openmp()


  if (iter_order > 0 .and. won==2) then
     write(*,'(a,i3,a)') 'WARNING: iter_order>0 (currently ',iter_order,')'
     write(*,*)          '         can not be used with won=2 (pixel-based weights)'
  endif

  !-----------------------------------------------------------------------
  !              allocate space for arrays
  !-----------------------------------------------------------------------
  polar = 0
  if (polarisation) polar = 1
  n_pols = 1 + 2*polar ! either 1 or 3

  ALLOCATE(map_TQU(0:npixtot-1,1:n_pols),stat = status)
  call assert_alloc(status,code,"map_TQU")

  if (twomaps) then
     ALLOCATE(map_TQU2(0:npixtot-1,1:n_pols),stat = status)
     call assert_alloc(status,code,"map_TQU2")
  endif

  ALLOCATE(w8ring_TQU(1:2*nsmax,1:n_pols),stat = status)
  call assert_alloc(status,code,"w8ring_TQU")

  ALLOCATE(dw8(1:2*nsmax,1:n_pols),stat = status)
  call assert_alloc(status,code,"dw8")

  ALLOCATE(fact_norm(0:nlmax),stat = status)
  call assert_alloc(status,code,"fact_norm")

  if (n_plm/=0) then
     ALLOCATE(plm(0:n_plm-1,1:plm_pol),stat = status)
     call assert_alloc(status,code,"plm")
  endif
  !-----------------------------------------------------------------------
  !                      reads the map
  !-----------------------------------------------------------------------
  PRINT *,"      "//code//"> Inputting map "

  fmissval = 0.0
  if (lowlreg > 0) fmissval = real(HPX_DBADVAL,kind=KMAP)
  call input_map(infile(1), map_TQU(0:,1:n_pols), &
       &  npixtot, n_pols, fmissval=fmissval, units=units_map(1:n_pols,1))
  units_map(1,1) = adjustl(units_map(1,1))
  units_pow(1) = trim(units_map(1,1))//'^2'

  if (twomaps) then
     call input_map(infile(2), map_TQU2(0:,1:n_pols), &
          &  npixtot, n_pols, fmissval=fmissval, units=units_map(1:n_pols,2))
     units_map(1,2) = adjustl(units_map(1,2))
     if (units_map(1,2) /= units_map(1,2)) then
        PRINT *,"      "//code//"> Warning: the 2 maps have different units"
     endif
  endif


  !-----------------------------------------------------------------------
  !                 reorder maps to RING if necessary
  !-----------------------------------------------------------------------
  if (order_map(1) == 2 .or. (twomaps .and. order_map(2) == 2)) then
     PRINT *,"      "//code//"> Convert Nest -> Ring "
     if (order_map(1) == 2) then
        call convert_nest2ring (nsmax, map_TQU)
        order_map(1) = 1
     endif
     if (twomaps .and. order_map(2) == 2) then
        call convert_nest2ring (nsmax, map_TQU2)
        order_map(2) = 1
     endif
  endif

  !-----------------------------------------------------------------------
  !                      reads the mask
  !-----------------------------------------------------------------------
  if (do_mask) then
     fmiss_mask = 0.0_KMAP ! set to mask=0 missing pixels
     pessim_mask = .true.
     PRINT *,"      "//code//"> Inputting mask from "//trim(maskfile)
     npix_mask = getsize_fits(maskfile, ordering=order_mask, nside=nside_mask, nmaps=nmasks)

     if (nside_mask <= 0) then
        print*,"Keyword NSIDE not found in FITS header!"
        call fatal_error(code)
     endif
     if ((order_mask /= 1).and.(order_mask /= 2)) then
        print*,"The ordering scheme of the mask must be RING or NESTED."
        print*,"No ordering specification is given in the FITS header!"
        call fatal_error(code)
     endif

     if (nside_mask == nsmax) then
        ! mask has same resolution as map: read and proceed
        allocate(mask(0:npix_mask-1,1:nmasks),stat=status)
        call assert_alloc(status,code,"mask")
        call input_map(maskfile, mask, npix_mask, nmasks, fmissval = fmiss_mask)
     else
        ! mask has different resolution than map: read
        allocate(mask(0:npixtot-1,1:nmasks),stat=status)
        call assert_alloc(status,code,"mask")
        allocate(mask_tmp(0:npix_mask-1,1:nmasks),stat=status)
        call assert_alloc(status,code,"mask_tmp")
        call input_map(maskfile, mask_tmp, npix_mask, nmasks, fmissval = fmiss_mask)
        ! then degrade or prograde
        print*,'Modifying mask to match map resolution:'
        write(*,'(a,i6,a,i6)') 'Nside:  ',nside_mask,' --> ',nsmax
        do i_mask = 1, nmasks
           pmask_tmp => mask_tmp(:,i_mask:i_mask)
           pmask     => mask(:,i_mask:i_mask)
           if (order_mask == 2) then ! nested scheme
              call udgrade_nest(pmask_tmp, nside_mask, pmask, nsmax, &
                   &            fmissval= fmiss_mask, pessimistic=pessim_mask)
           else ! ring scheme
              call udgrade_ring(pmask_tmp, nside_mask, pmask, nsmax, &
                   &            fmissval= fmiss_mask, pessimistic=pessim_mask)
           endif
        enddo
        deallocate(mask_tmp)
        nside_mask = nsmax ; npix_mask = npixtot
     endif

     !--------------------------------------
     ! **** what if nmasks /= nmaps ??? ****
     !---------------------------------------
     do i_mask = 1, nmasks
        max_mask = maxval(mask(:,i_mask))
        min_mask = minval(mask(:,i_mask))
        if (max_mask > 0. .and. min_mask >=0.) then
           write(*,'(a,i2,a)') 'Field #', i_mask,' is OK'
           exit
        endif
        write(*,'(a,i2,a,2(1pg11.3))') 'Field #', i_mask, &
             & ' of '//trim(maskfile)//' has range ', min_mask, max_mask
        print*,' not a valid mask, will try next field'
        if (i_mask == nmasks) then
           print*,'Error : No valid mask found in '//trim(maskfile)
           call fatal_error(code)
        endif
     enddo

     ! keep only first valid mask found
     pmask => mask(:,i_mask:i_mask)
     nmasks = 1

     ! make sure mask and sky map have same ordering
     if (order_mask /= order_map(1)) then 
        PRINT *,"      "//code//"> Reorder mask "      
        if (order_mask == 2) call convert_nest2ring (nside_mask, pmask)
!!!!        if (order_mask == 1) call convert_ring2nest (nside_mask, pmask)
        order_mask = order_map(1)
     endif

  else ! create mask place holder
     allocate(pmask(0:1,1))
  endif 

  !-----------------------------------------------------------------------
  !                      remove dipole
  !-----------------------------------------------------------------------

  if (lowlreg > 0) then
     PRINT *,"      "//code//"> Remove monopole (and dipole) from Temperature map"
     call remove_dipole(nsmax, map_TQU(0:npixtot-1,1), order_map(1), lowlreg, &
          & mono_dip(0:,1), zbounds, fmissval=fmissval, mask=pmask(:,1))
     if (twomaps) then
        call remove_dipole(nsmax, map_TQU2(0:npixtot-1,1), order_map(2), lowlreg, &
             & mono_dip(0:,2), zbounds, fmissval=fmissval, mask=pmask(:,1))
     endif

     if (fmissval /= 0.0) then
        do j=1,n_pols
           do ipix=0, npixtot-1
              if (abs(map_TQU(ipix,j)-fmissval) <= abs(1.e-6*fmissval)) map_TQU(ipix,j) = 0.0_kmap
           enddo
           if (twomaps) then
              do ipix=0, npixtot-1
                 if (abs(map_TQU2(ipix,j)-fmissval) <= abs(1.e-6*fmissval)) map_TQU2(ipix,j) = 0.0_kmap
              enddo
           endif
        enddo
     endif

     do i=1, n_maps
        write(unit=*,fmt=8000)  " Monopole = ", mono_dip(0,i)," "//trim(units_map(1,i))
        if (lowlreg > 1) then
           call vec2ang( mono_dip(1:3,i), theta_dip, phi_dip)
           ampl_dip = sqrt(sum(mono_dip(1:3,i)**2))
           long_dip =      phi_dip  /PI*180.0_dp
           lat_dip  = 90.0-theta_dip/PI*180.0_dp
           write(unit=*,fmt=8000)  " Dipole   = ",ampl_dip," "//trim(units_map(1,i))
           write(unit=*,fmt=8002) "(long.,lat.) = (",long_dip,lat_dip,") Deg"
        endif
     enddo
8000 format(a,1pg13.3,a)
8002 format(a,1pg10.3,', ',1pg10.3,a)

  endif
  !-----------------------------------------------------------------------
  !                      Input weights
  !-----------------------------------------------------------------------

  dw8=0.0_dp ! read as DP
  if (trim(infile_w8)/="") then
     if (won == 1) then
        ! ring based weights: treated as such in map2alm
        n_rings = 2 * nsmax
        n_rings_read = getsize_fits(infile_w8, nmaps=nmw8)

        if (n_rings_read /= n_rings) then
           PRINT *," "
           print *,"wrong ring weight file:"//trim(infile_w8)
           call fatal_error(code)
        endif
        nmw8 = min(nmw8, n_pols)
        
        PRINT *,"      "//code//"> Inputting "//trim(strw8file)
        call input_map(infile_w8, dw8, n_rings, nmw8, fmissval=0.0_dp)
     endif

     if (won == 2) then 
        ! pixel based weights: turned into full sky map and multiply with mask
        allocate(pw8map (0:npixtot-1,   1:1), stat=status)
        call assert_alloc(status,code,"pw8map")

        PRINT *,"      "//code//"> Inputting and unfolding "//trim(strw8file)

!         n_pw8_read = getsize_fits(infile_w8)
!         allocate(pw8list(0:n_pw8_read-1,1:1), stat=status)
!         call assert_alloc(status,code,"pw8list")
!         call input_map(infile_w8, pw8list, n_pw8_read, 1)
!         call unfold_weights(nsmax, won, pw8list(:,1), pw8map(:,1))
!         deallocate(pw8list)

        call unfold_weightsfile(infile_w8, pw8map(:,1))

        if (do_mask) then
           do i=1,size(pmask,2)
              pmask(1:npixtot,i:i) = pmask(1:npixtot,i:i) * pw8map(0:npixtot-1,1:1)
           enddo
        else
           deallocate(pmask)
           pmask => pw8map(:,1:1)
        endif
        !print*,pmask(1:12,1)
        !print*,pmask(npixtot-11:npixtot,1)
     endif
  endif
  w8ring_TQU = 1.0_dp + dw8

  !-----------------------------------------------------------------------
  !                    precomputed plms
  !-----------------------------------------------------------------------

  if (n_plm/=0) then
     if (plm_pol.eq.1) then
        PRINT*,"      "//code//"> Reading precomputed scalar P_lms "
     else
        PRINT*,"      "//code//"> Reading precomputed tensor P_lms "
     endif
     !call read_dbintab(infile_plm,plm,n_plm,plm_pol,nullval,anynull=bad) ! replaced 2010-11-24
     call read_bintab(infile_plm, plm, n_plm, plm_pol, nullval, anynull=bad)
     if (bad) call fatal_error("Missing points in P_lm-file!")
  endif

  !-----------------------------------------------------------------------
  !                        map to alm
  !  map2alm_iterative will apply mask to map
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

  if (twomaps) then
     allocate(alm_TGC2(1:n_pols, 0:nlmax, 0:nmmax),stat = status)
     call assert_alloc(status,code,"alm_TGC2")
     if (n_plm/=0) then
        call map2alm_iterative(nsmax, nlmax, nmmax, iter_order, map_TQU2, alm_TGC2, &
          &                 zbounds=zbounds, w8ring=w8ring_TQU, plm=plm, mask=pmask)
     else
        call map2alm_iterative(nsmax, nlmax, nmmax, iter_order, map_TQU2, alm_TGC2, &
          &                 zbounds=zbounds, w8ring=w8ring_TQU, mask=pmask)
     endif
     deallocate( map_TQU2 )
  endif

  deallocate( map_TQU )
  deallocate(w8ring_TQU)
  deallocate(dw8)
  if (allocated(pw8map)) deallocate(pw8map)

  !-----------------------------------------------------------------------
  !                  power spectrum calculation
  !-----------------------------------------------------------------------

  PRINT *,"      "//code//"> Compute power spectrum "

  ncl = 1
  if (polarisation) ncl = 4
  if (polarisation .and. bcoupling)    ncl = 6
  asym_cl = polarisation .and. bcoupling .and. twomaps
  if (asym_cl) ncl = 9
  allocate(clout(0:nlmax,1:ncl))

  if (twomaps) then
     call alm2cl(nlmax, nmmax, alm_TGC, alm_TGC2, clout)
  else
     call alm2cl(nlmax, nmmax, alm_TGC, alm_TGC, clout)
  endif
  !-----------------------------------------------------------------------
  !                  outputs the power spectrum
  !-----------------------------------------------------------------------

  if (trim(outfile)/="") then
     PRINT *,"      "//code//"> outputting Power Spectrum "

!      do l = 0, nlmax
!         fact_norm(l) = 1.0_KCL !  Normalisation=1, Output 'pure' C_ls
!      enddo
!      do i = 1, ncl
!         clout(:,i) = clout(:,i) * fact_norm(:)
!      enddo

!     units_pow(1) = "'"//trim(units_pow(1))//"'" ! add quotes to keep words together
     call write_minimal_header(header,'cl', &
          creator = CODE, version = VERSION, &
          polar = polarisation, bcross = bcoupling, asym_cl = asym_cl, &
          nlmax = nlmax, nside = nsmax, coordsys = coordsys, units = units_pow(1) )
     if (twomaps) then
        call add_card(header,"EXTNAME","'ANALYSED CROSS POWER SPECTRUM'",update=.true.)
     else
        call add_card(header,"EXTNAME","'ANALYSED AUTO POWER SPECTRUM'",update=.true.)
     endif
     call add_card(header,"CUT-SKY", theta_cut_deg," [deg] Symmetric Cut Sky")
     call add_card(header,"NITERALM",iter_order,   " Number of iterations for a_lms extraction")
     call add_card(header,"LREGRESS",lowlreg,      " Regression of low multipoles (0/1/2)")
     call add_card(header,"NMAPS_IN",n_maps,       " Number of maps (co-)analysed (1/2)")
     call add_card(header,"QUADSCHM",won,          " Quadrature scheme (0/1/2)")
     call add_card(header)
     if (twomaps) then
        call add_card(header,"HISTORY","Input Map1 = "//TRIM(infile(1)))
        call add_card(header,"HISTORY","Input Map2 = "//TRIM(infile(2)))        
     else
        call add_card(header,"HISTORY","Input Map = "//TRIM(infile(1)))
     endif
     nlheader = SIZE(header)
     !print*,'wrtting file'
     call write_asctab (clout, nlmax, ncl, header, nlheader, outfile)
  endif
  !-----------------------------------------------------------------------
  !                  outputs the alm
  !-----------------------------------------------------------------------
  do j = 1,n_maps
     if (trim(outfile_alms(j))/="") then

        PRINT *,"      "//code//"> outputting alms "
!        units_map(1,j) = ''''//trim(units_map(1,j))//'''' ! add quotes to keep words together

        ! write each component (T,G,C) in a different extension of the same file
        do i=1,1+polar*2
           call write_minimal_header(header,'alm', &
                creator = CODE, version = VERSION, &
                nside = nsmax, coordsys = coordsys, &
                nlmax = nlmax, nmmax = nmmax, polar = polarisation, &
                units = units_map(1,j)) 
           if (i == 1) then
              call add_card(header,"EXTNAME","'ANALYSED a_lms (TEMPERATURE)'")
           elseif (i == 2) then
              call add_card(header,"EXTNAME","'ANALYSED a_lms (GRAD / ELECTRIC component)'")
           elseif (i == 3) then
              call add_card(header,"EXTNAME","'ANALYSED a_lms (CURL / MAGNETIC component)'")
           endif
           call add_card(header,"CUT-SKY", theta_cut_deg," [deg] Symmetric Cut Sky")
           call add_card(header,"NITERALM",iter_order,   " Number of iterations for a_lms extraction")
           call add_card(header,"LREGRESS",lowlreg,      " Regression of low multipoles (0/1/2)")
           call add_card(header,"HISTORY","Input Map = "//TRIM(infile(j)))
           call add_card(header,"QUADSCHM",won,          " Quadrature scheme (0/1/2)")
           nlheader = SIZE(header)
           if (j == 1) call dump_alms(outfile_alms(j), alm_TGC (i,0:nlmax,0:nmmax), &
                &                     nlmax, header, nlheader, i-1_i4b)
           if (j == 2) call dump_alms(outfile_alms(j), alm_TGC2(i,0:nlmax,0:nmmax), &
                &                     nlmax, header, nlheader, i-1_i4b)
        enddo ! i: extension
     endif
  enddo ! j: file

  !-----------------------------------------------------------------------
  !                      deallocate memory for arrays
  !-----------------------------------------------------------------------

  deallocate( alm_TGC )
  if (twomaps) deallocate( alm_TGC2 )

  !-----------------------------------------------------------------------
  !                      report card
  !-----------------------------------------------------------------------
  call wall_clock_time(time1)
  call cpu_time(ptime1)
  clock_time = time1 - time0
  ptime      = ptime1 - ptime0

  write(*,9000) " "
  write(*,9000) " Report Card for "//code//" analysis run"
  write(*,9000) "----------------------------------------"
  write(*,9000) " "
  if (.not. polarisation) then
     write(*,9000) "      Temperature alone"
  else
     write(*,9000) "    Temperature + Polarisation"
  endif
  write(*,9000) " "
  write(*,9000) " Input map              : "//TRIM(infile(1))
  if (trim(infile(2))/="") write(*,9000) " Input map 2            : "//TRIM(infile(2))
  write(*,9006) " Multipole range        : ",lowlreg ," <= l <= ", nlmax
  write(*,9010) " Number of pixels       : ", npixtot
  write(*,9020) " Pixel size in arcmin   : ", pix_size_arcmin
  write(*,9020) " Symmetric cut in deg   : ", theta_cut_deg
  if (do_mask) write(*,9000) &
       &        " Mask read from         : "//trim(maskfile)
  write(*,9000) " Output power spectrum  : "//TRIM(outfile)
  if (outfile_alms(1)/="") write(*,9000) " Output alm             : "//TRIM(outfile_alms(1))
  if (outfile_alms(2)/="") write(*,9000) " Output alm 2           : "//TRIM(outfile_alms(2))
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

