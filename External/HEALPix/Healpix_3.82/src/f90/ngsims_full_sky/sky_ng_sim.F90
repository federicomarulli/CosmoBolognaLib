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

!******************************************************
! Sarah Smith (sjm84@mrao.cam.ac.uk)
! and Graca Rocha (graca@mrao.cam.ac.uk)
! June 2004
! Module sub_ngpdf_sho based on code by Michael
! Hobson and Graca Rocha
! Module sub_ngpdf_powergauss based on code by
! Anthony Challinor
! Please include an appropriate acknowledgement 
! in any publications based on work that has
! made use of this package - 'NGsims'
!******************************************************
Program sky_ng_sim
  ! Makes a full-sky non-Gaussian simulation using
  ! values drawn from pdf made from eigenstates of
  ! a simple harmonic oscillator or powers of a Gaussian
  ! Uses Healpix pixelisation and subroutines
  ! Some of the program is based on the Healpix synfast program

  USE healpix_types
  USE alm_tools, ONLY : map2alm_iterative, alm2map, alm2map_der, pow2alm_units
  USE fitstools, ONLY : read_asctab, write_bintab, dump_alms
  USE pix_tools, ONLY : nside2npix
  USE head_fits, ONLY : add_card, merge_headers, get_card, write_minimal_header, del_card
!  USE utilities, ONLY : die_alloc
  use misc_utils, only : wall_clock_time, assert_alloc, brag_openmp, fatal_error
  USE extension, ONLY : getArgument, nArguments
  USE paramfile_io, ONLY : paramfile_handle, parse_int, parse_init, &
       & parse_string, parse_lgt, parse_double, parse_check_unused, &
       & concatnl, scan_directories, parse_summarize, parse_finish, &
       & get_healpix_data_dir, get_healpix_test_dir,get_healpix_pixel_window_file 

  USE sub_ngpdf_sho
  USE sub_ngpdf_powergauss
  USE sky_sub

  Implicit none

  INTEGER(I4B) :: nsmax ! The value of N_{side} in the pixelisation scheme
  INTEGER(I4B) :: nlmax, lmax ! The maximum l-value used
  INTEGER(I4B) :: nmmax ! The maximum m-value used (set equal to nlmax)

  COMPLEX(SPC), DIMENSION(:,:,:), ALLOCATABLE :: alm_T !The a_{lm} values
  REAL(SP),     DIMENSION(:, :),  ALLOCATABLE :: map_T !The pixel values in real space
  REAL(SP), DIMENSION(:,:), ALLOCATABLE :: cl_T !The Cl values
  REAL(KIND=DP),     DIMENSION(:,:),   ALLOCATABLE :: w8ring_T
  REAL(SP),     DIMENSION(:,:),     ALLOCATABLE :: tmp_2d

!  INTEGER(I4B), DIMENSION(8,2) :: values_time
!  REAL(SP) :: clock_time
  real(SP) :: time0, time1, clock_time, ptime0, ptime1, ptime
  INTEGER(I4B) :: status

  INTEGER(I4B) npixtot !The total number of pixels
  CHARACTER(LEN=filenamelen)          :: parafile = ''
  CHARACTER(LEN=filenamelen)          :: infile
  character(len=FILENAMELEN)          :: outfile_alms
  CHARACTER(LEN=filenamelen)          :: outfile
  CHARACTER(LEN=filenamelen)          :: windowfile
  CHARACTER(LEN=filenamelen)          :: windowname
  CHARACTER(LEN=filenamelen)          :: def_dir, def_file
  CHARACTER(LEN=filenamelen)          :: usr_dir, usr_file
  CHARACTER(LEN=filenamelen)          :: final_file
  CHARACTER(LEN=filenamelen)          :: healpixtestdir
  Character(LEN=filenamelen)          :: beam_file
  character(len=filenamelen)          :: description
  character(len=100)                  :: chline
  LOGICAL(LGT) :: ok, fitscl, polarisation = .FALSE.
  logical(LGT) :: do_map, output_alms
  CHARACTER(LEN=80), DIMENSION(1:180) :: header, header_PS !, header_file
  CHARACTER(LEN=*), PARAMETER :: code = "sky_ng_sim"
  character(len=*), parameter :: VERSION = HEALPIX_VERSION
  character(len=80), dimension(1:1) :: units_power, units_map
  CHARACTER(LEN=20)                            ::  string_quad
!  CHARACTER(LEN=80) :: temptype, com_tt
!  REAL(SP) ::  quadrupole
!  INTEGER(I4B) :: count_tt
  INTEGER(I4B) nlheader
  Integer :: i,m, n_pols, n_maps, simul_type, deriv
  Real(DP) :: sigma0, factor
  Real(DP),dimension(1:3) :: nu !Added to match with f90 version of shodev_driver
  !nu(i) is the ith moment of the distribution. 
  !The size of nu is fixed to 3 to prevent problems with passing unallocated arrays
  Real(DP) :: power, rms_alm, mean
  Integer :: count, iter_order
  Integer(I4B) :: pdf_choice !Used to chose type of pdf required
  logical :: plot
  !  External shodev_driver
  Real(SP) :: Tmin, Tmax, convert
  Integer, parameter :: npts = 1000
  Real(SP) :: xline(npts), yline(npts), xstep, xval
  type(paramfile_handle) :: handle
  Real(DP) :: fwhm_arcmin
  Real(SP) :: fwhm_deg
  Character(len=30) :: xlabel
  !-----------------------------------------------------------------------
  !                    get input parameters and create arrays
  !-----------------------------------------------------------------------

  call wall_clock_time(time0)
  call cpu_time(ptime0)
  !     --- read parameters interactively if no command-line arguments
  !     --- are given, otherwise interpret first command-line argument
  !     --- as parameter file name and read parameters from it:
  if (nArguments() == 0) then
     parafile=''
  else
     if (nArguments() /= 1) then
        print *, " Usage: "//code//" [parameter file name]"
        call fatal_error()
     end if
     call getArgument(1,parafile)
  end if
  handle = parse_init(parafile)

  PRINT *, " "
  PRINT *,"                    "//code//" "//version
  Write(*, '(a)')  "*** Simulation of a non-Gaussian full-sky temperature map ***"

  !     --- choose temp. only or  temp. + deriv. ---
  description = concatnl( &
       & " Do you want to simulate", &
       & " 1) Temperature only", &
       !& " 2) Temperature + polarisation", &
       & " 3) Temperature + 1st derivatives", &
       & " 4) Temperature + 1st & 2nd derivatives")
  simul_type = parse_int(handle, 'simul_type', default=1, valid=(/1,3,4/), descr=description)
  deriv = 0 ! no derivatives
  if (simul_type == 3) deriv = 1
  if (simul_type == 4) deriv = 2
  n_pols = 1
  n_maps = max(1, 3*deriv)

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


  !     --- get filename for input power spectrum ---
  chline = ''
  healpixtestdir = get_healpix_test_dir()
  if (trim(healpixtestdir)/='') chline = trim(healpixtestdir)//'/cl.fits'
  description = concatnl( "", &
       & " Enter input Power Spectrum filename", &
       & " Can be either in FITS format or in the .dat format produced by CAMB.")
  infile = parse_string(handle, 'infile', default=chline, descr=description, filestatus='old')
  ! Find out whether input cl file is a fits file
  fitscl = (index(infile, '.fits') /= 0) 

  !     --- gets the output sky map filename ---
  description = concatnl(&
       & "", &
       & " Enter Output map file name (eg, test.fits) :", &
       & "  (or !test.fits to overwrite an existing file)" )
  outfile = parse_string(handle, "outfile", &
       default="!test.fits", descr=description, filestatus="new")
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

  !     --- gets the fwhm of beam ---
  description = concatnl(&
       & "", &
       & " Enter FWMH of beam in arcminutes (0 for no beam):")
  fwhm_arcmin = parse_double(handle, "fwhm_arcmin", &
       default=0d0, descr=description, vmin = 0d0)

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
     print*,'The beam file '//trim(beam_file)//' will be used instead.'
  endif

  ! including pixel window function, EH-2008-03-05
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
     usr_dir = parse_string(handle,'winfiledir',default='',descr=description)
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

  PRINT *," "

  !-----------------------------------------------------------------------

  nmmax   = nlmax
  npixtot = nside2npix(nsmax)

  !-----------------------------------------------------------------------
  !                  allocates space for arrays
  !-----------------------------------------------------------------------

  !  ALLOCATE(units_alm(1:1),units_map(1:1+2*polar),stat = status)
  !  call assert_alloc(status, code,"units_alm & units_map") 

  ALLOCATE(alm_T(1:1,0:nlmax, 0:nmmax),stat = status)
  call assert_alloc(status, code,"alm_T")

  ALLOCATE(map_T(0:npixtot-1, 1:n_maps),stat = status)
  call assert_alloc(status, code,"map_T")

  ALLOCATE(w8ring_T(1:2*nsmax,1:1),stat = status)
  call assert_alloc(status, code,"w8ring_T")

  ! For now, not using ring weights for quadrature correction
  w8ring_T = 1.d0

  ! single analysis
  iter_order = 0 

  ALLOCATE(cl_T(0:nlmax,1:1),stat = status)
  call assert_alloc(status, code,"cl_T")

  !------------------------------------------------------------------------
  ! Read in the input power spectrum
  !-----------------------------------------------------------------------

  cl_T = 0.0
  !New lines added 8th June 2004
  lmax = nlmax
  call read_powerspec(infile, nsmax, lmax, cl_T, header_PS, fwhm_arcmin, units_power, beam_file = beam_file, winfile = windowfile) 
  call pow2alm_units(units_power, units_map)
  call del_card(header_PS, (/ "TUNIT#","TTYPE#"/)) ! remove TUNIT* and TTYPE* from header to avoid confusion later on

  !------------------------------------------------------------------------
  ! Draw pixel values in real space from non-Gaussian distribution
  !------------------------------------------------------------------------

  Write (*,*) "Creating non-Gaussian map with flat power spectrum"

  !     --- Chose the type of pdf to use ---
  description = concatnl(&
       & "", &
       & " Select non-Gaussian pdf to use: Simple harmonic oscillator (1)", &
       & "or powers of a Gaussian (2)" )
  pdf_choice = parse_int(handle, 'pdf_choice', default=1, vmin = 1, vmax = 2, descr=description)  

  !Write (*,*) "Use SHO pdf (0) or powers of a Gaussian pdf (1) ?"
  !Read (*,*) m
  If (pdf_choice .eq. 1) Then
     call shodev_driver(handle, npixtot, sigma0, map_T(:,1), nu)
  Else
     !     call powergauss_driver(npixtot, npixtot, sigma0, map_T, nu)
     call powergauss_driver(handle, npixtot, sigma0, map_T(:,1), nu)
  End If
  print*,'NG map MIN and MAX:',minval(map_T(:,1)),maxval(map_T(:,1))
  !Normalise the map
  If (nu(1) .ge. 0) Then !nu (1) is set to -1 if highest non-zero alpha > 3
     ! ie if we can't calculate moments analytically
     map_T = map_T / sqrt(nu(2))

     !Compute the rms value of the pixels
     power = 0
     Do i = 0, npixtot-1
        power = power + map_T(i,1)**2
     End Do
     power = Sqrt(power/npixtot)
     Write (*,*) "rms value of the pixels is ", power
  Else
     Write (*,*) "Can't calculate theoretical value of nu(2) to normalise"
     Write (*,*) "Normalising using rms pixel value instead"
     Write (*,*) "note - this method may change the statistical properties of the map"
     !Compute the mean value of the pixels
     mean = SUM(map_T(:,1)*1.d0)/npixtot
     Write (*,*) 'Mean pixel value is ',mean,' - adjusting to zero'
     map_T(:,1) = map_T(:,1) - mean
     !Compute the rms value of the pixels
     power = 0
     Do i = 0, npixtot-1
        power = power + map_T(i,1)**2
     End Do
     power = Sqrt(power/npixtot)
     Write (*,*) "rms value of the pixels is ", power ," - setting to 1.0"
     map_T(:,1) = map_T(:,1) / power    
     !Check
     power = 0
     Do i = 0, npixtot-1
        power = power + map_T(i,1)**2
     End Do
     power = Sqrt(power/npixtot)
     Write (*,*) "Now rms value of the pixels is ", power
  End If
  !------------------------------------------------------------------------
  ! Now compute the a_{lm} values from this distribution
  !------------------------------------------------------------------------

  Write (*,*) "Computing the values of the a_{lm}"
!  call map2alm(nsmax, nlmax, nmmax, map_T, alm_T, -1000.d0, w8ring_T)
!  call map2alm(nsmax, nlmax, nmmax, map_T, alm_T, (/ 0.d0, 0.d0 /), w8ring_T)
  call map2alm_iterative(nsmax, nlmax, nmmax, iter_order, map_T(:,1:1), alm_T, w8ring=w8ring_T)


  !Compute the rms value of the a_{lm}
  power = 0.0
  count = 0
  Do i = 1, nlmax
     power = power + abs(alm_T(1,i,0))**2.
     count = count + 1
     Do m = 1, i
        power = power + 2.0*abs(alm_T(1,i,m))**2.
        count = count + 2
     End Do
  End Do
  rms_alm = sqrt(power/count)
  Write (*,*) "rms value of the a_{lm} is ",rms_alm
  factor = sqrt(FOURPI/npixtot)
  Write (*,*) "Expected rms value is ", factor

  !------------------------------------------------------------------------
  ! Multiply the a_{lm} by the correct power spectrum
  !------------------------------------------------------------------------

  Do i = 1, nlmax
     !Write (*,*) 'cl(',i,') is ', cl_T(i,1)
     alm_T(1,i,0:i) = alm_T(1,i,0:i)*sqrt(cl_T(i,1))/factor
  End Do

  !------------------------------------------------------------------------
  ! Now transform back to a real space map
  !------------------------------------------------------------------------

  select case (deriv)
  case(0) ! temperature only
     call alm2map(    nsmax,nlmax,nmmax,alm_T,map_T(:,1))
  case(1)
     call alm2map_der(nsmax,nlmax,nmmax,alm_T,map_T(:,1),map_T(:,2:3))
  case(2)
     call alm2map_der(nsmax,nlmax,nmmax,alm_T,map_T(:,1),map_T(:,2:3),map_T(:,4:6))
  case default
     print*,'Not valid case'
     print*,'Derivatives: ', deriv
     call fatal_error(code)
  end select

  !------------------------------------------------------------------------
  ! Plot histogram of pixel distribution (if required)
  !------------------------------------------------------------------------
  
#ifdef PGPLOT
  description = concatnl( "", &
       & " Plot histogram of pixel values? (True or False)")
  plot = parse_lgt(handle, 'plot', default=.false., descr=description)
!  Write (*,*) 'Plot histrogram of pixel values (0) or not (1)?'
!  Read (*,*) m
  If (plot) Then
     Tmin = 0.0
     Tmax = 0.0
     do i=1,npixtot
        if (map_T(i,1) .lt. Tmin) Tmin=map_T(i,1)
        if (map_T(i,1) .gt. Tmax) Tmax=map_T(i,1)
     end do
     call pgbegin(0,'?',1,1)
     if (.not. fitscl) then
        ! If reading in from .dat file Cls are already converted to uK
        convert = 1.0
        xlabel = 'Temperature / \gmK'
     else if (ABS(cl_T(2,1) - 1.0) .lt. 0.05) then
        ! If values in fits file are normalised to Cl(2) = 1.0
        ! convert to uK for plot, assuming COBE normalisation
        ! of Qrms = 18uK
        convert = sqrt(4.0*PI/5.0)*18.0
        xlabel = 'Temperature / \gmK'
     else
        ! Otherwise, units unknown
        convert = 1.0
        xlabel = 'Pixel value'
     end if
     Tmin = Tmin * convert
     Tmax = Tmax * convert
     call pghist(npixtot,map_T*convert,Tmin,Tmax,200,0)
     !Pixel values scaled to convert to uK if input fits file of cl values
     !which are normalised to C_2 = 1.0
     call pglab(xlabel,'Number of pixels','')
     !Calculate 2*mean square value of T (ie 2 sigma^2)
     power = 0
     Do i = 1, npixtot
        power = power + map_T(i,1)**2
     End Do
     power = 2 * power*(convert**2) / npixtot
     xstep=(Tmax-Tmin)/real(npts-1)
     do i=1,npts
        xval=Tmin+real(i)*xstep
        xline(i)=xval
        yline(i)=exp(-xval**2/power)/sqrt(PI*power)
        yline(i)=yline(i)*real(npixtot)*(Tmax-Tmin)/200.
     end do
     call pgsci(13)
     call pgline(npts, xline, yline)
     call pgend
  End If   
#endif

  call parse_check_unused(handle, code=code)
  call parse_summarize(handle, code=code)
  call parse_finish(handle)
  call brag_openmp()

  !-----------------------------------------------------------------------
  !                      write the alm to FITS file
  !-----------------------------------------------------------------------
  if (output_alms) then

     PRINT *,"      "//code//"> outputting alms "

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
             nlmax = nlmax, nmmax = nmmax, &
             !randseed = ioriginseed, &
             units = units_map(i), polar = polarisation)
        if (i == 1) then
           call add_card(header,"EXTNAME","'SIMULATED a_lms (TEMPERATURE)'", update = .true.)
        elseif (i == 2) then
           call add_card(header,"EXTNAME","'SIMULATED a_lms (GRAD / ELECTRIC component)'", update = .true.)
        elseif (i == 3) then
           call add_card(header,"EXTNAME","'SIMULATED a_lms (CURL / MAGNETIC component)'", update = .true.)
        endif
        !if (input_cl) then
           call add_card(header,"HISTORY","File with input C(l):")
           call add_card(header,"HISTORY",trim(infile))
           call add_card(header,"HISTORY","These alms are multiplied by pixel and beam window functions")
           call add_card(header,"NSIDE"   ,nsmax,   "Resolution parameter for HEALPIX")
           if (trim(beam_file) == '') then
              call add_card(header,"FWHM"    ,fwhm_deg   ," [deg] FWHM of gaussian symmetric beam")
           else
              call add_card(header,"BEAM_LEG",trim(beam_file), &
                   & "File containing Legendre transform of symmetric beam")
           endif
           call add_card(header,"PDF_TYPE",pdf_choice,"1: Harmon. Oscill. ;2: Power of Gauss.")
        !endif
        call add_card(header,"HISTORY")
        nlheader = SIZE(header)
        call dump_alms (outfile_alms,alm_T(i,0:nlmax,0:nmmax),nlmax,header,nlheader,i-1_i4b)
     enddo
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
       nside = nsmax, ordering = 'RING', &   !, coordsys = coordsys, &
       fwhm_degree = fwhm_arcmin / 60.d0, &
       beam_leg = trim(beam_file), &
       polar = polarisation, &
       deriv = deriv, &
       creator = CODE, version = VERSION, &
       nlmax = lmax, &
       !randseed = ioriginseed, &
       units = units_map(1) )
  !Note using lmax rather than nlmax - lmax will be smaller than nlmax if not enough values
  ! in input power spectrum file

  call add_card(header,"PDF_TYPE",pdf_choice,"1: Harmon. Oscill. ;2: Power of Gauss.")
  call add_card(header,"EXTNAME","'SIMULATED MAP'", update=.true.)
  call add_card(header,"COMMENT","*************************************")

!   allocate(tmp_2d(0:npixtot-1,1:1))
!   tmp_2d(:,1) = map_T
!   call write_bintab(tmp_2d, npixtot, 1, header, nlheader, outfile)
!   deallocate(tmp_2d)
  call write_bintab(map_T, npixtot, n_maps, header, nlheader, outfile)
!!$  endif

  !-----------------------------------------------------------------------
  !                      deallocate memory for arrays
  !-----------------------------------------------------------------------

  DEALLOCATE( map_T )
  DEALLOCATE( alm_T )
  !  DEALLOCATE( units_alm, units_map )

  !-----------------------------------------------------------------------
  !                      output and report card
  !-----------------------------------------------------------------------
  call wall_clock_time(time1)
  call cpu_time(ptime1)
  clock_time = time1 - time0
  ptime      = ptime1 - ptime0

  WRITE(*,9000) " "
  WRITE(*,9000) " Report Card for "//code//" simulation run"
  WRITE(*,9000) "----------------------------------------"
  WRITE(*,9000) " "
  WRITE(*,9000) " Input power spectrum : "//TRIM(infile)
  WRITE(*,9010) " Multipole range      : 0 < l <= ", nlmax
  WRITE(*,9010) " Number of pixels     : ", npixtot
  !  WRITE(*,9020) " Pixel size in arcmin : ", pix_size_arcmin
  !  WRITE(*,9010) " Initial random # seed: ", ioriginseed
  if (trim(beam_file) == '') then
     write(*,9020) " Gauss. FWHM in arcmin: ", fwhm_arcmin
  else
     write(*,9000) " Beam file: "//trim(beam_file)
  endif
  WRITE(*,9000) " Output map           : "//TRIM(outfile)
  write(*,9030) " Clock and CPU time [s] : ", clock_time, ptime

  !-----------------------------------------------------------------------
  !                       end of routine
  !-----------------------------------------------------------------------

  WRITE(*,9000) " "
  WRITE(*,9000) " "//code//"> normal completion"

9000 format(a)
9010 format(a,i16)
9020 format(a,g20.5)
9030 format(a,f11.2,f11.2)


  Stop

End Program sky_ng_sim
