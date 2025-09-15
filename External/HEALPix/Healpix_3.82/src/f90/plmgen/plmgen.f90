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
!===================================================================================
program plmgen
  ! scalar implementation
  !===================================================================================
  ! precomputation of scalar and tensorial spherical harmonics
  ! version 1.0 : ?
  ! version 1.1 : ?
  ! version 1.2 :
  !   Aug 14, 2000, EH, Caltech, switch to CMBFAST normalisation of tensorial spherical harmonics
  !       (see Zaldarriaga, astro-ph/9709271; ApJ, 503, 1 (1998))
  ! version 1.3 : Aug 2004, EH, Ipac
  !      put calculations in plm_gen in alm_tools
  ! version 1.4: June 2008: check that job is compatible with 4 bytes variable,
  !              more verbose parsing
  !===================================================================================
  USE healpix_types
  USE misc_utils, ONLY: assert_alloc, wall_clock_time, strlowcase, brag_openmp, fatal_error
  USE head_fits, ONLY: add_card
  USE fitstools, ONLY: write_plm
  USE pix_tools, ONLY: nside2npix
  USE alm_tools, ONLY: plm_gen
  USE extension, ONLY : getArgument, nArguments
  USE paramfile_io, ONLY : paramfile_handle, parse_init, parse_int, &
         parse_string, concatnl, parse_check_unused, parse_summarize, parse_finish

  implicit none

  INTEGER(I4B) :: nsmax, nlmax, nmmax
  INTEGER(I4B) :: nd
  REAL(DP), DIMENSION(:,:),ALLOCATABLE :: plm
  LOGICAL(LGT) :: polarisation

  CHARACTER(LEN=filenamelen)        :: outfile
  CHARACTER(LEN=80), DIMENSION(160) :: header
  CHARACTER(LEN=*), PARAMETER :: code = 'PLMGEN'
  character(len=*), parameter :: VERSION = HEALPIX_VERSION

  INTEGER(I4B) :: n_plm
  integer(i8b) :: n_plm_8
  INTEGER(I4B) :: status, nlheader, simul_type

  type(paramfile_handle)        :: handle
  character(LEN=filenamelen)    :: parafile = ''
  character(len=filenamelen)    :: description
  character(len=100)            :: chline
  character(len=20)             :: lcode

  real(kind=SP)                 :: clock_time, time0, time1
  real(kind=SP)                 :: ptime, ptime0, ptime1

  !===================================================================================

  call wall_clock_time(time0)
  call cpu_time(ptime0)
  lcode = trim(strlowcase(code))

  if (nArguments() == 0) then
     parafile=''
  else
     if (nArguments() /= 1) then
        write(*,'(a)') " Usage:"//trim(lcode)//" [parameter file name]"
        stop 1
     end if
     call getArgument(1,parafile)
  end if

  PRINT *, " "
  PRINT *,"                    "//code//" "//version
  write(*,'(a)') &
       & " *** Precomputation of spherical harmonics (modulus of the complex Ylm) ***" 
  write(*,'(a)') &
       & " ***   used to analyse/synthetize Healpix maps of a given resolution    ***"
  PRINT *, " "

  handle = parse_init(parafile)

  !     --- choose temp. only or  temp. + pol. ---
  description = concatnl( &
       & " Do you want to compute the spherical harmonics for", &
       & " 1) Temperature only", &
       & " 2) Temperature + polarisation")
  simul_type = parse_int(handle, 'simul_type', default=1, vmin=1, vmax=2, descr=description)
  if (simul_type .eq. 1) polarisation = .false.
  if (simul_type .eq. 2) polarisation = .true.

  !     --- get resolution parameter ---
3 continue
  description = concatnl( &
       & "", &
       & " Enter the resolution parameter (Nside) of the map to be analysed/synthetized: ",&
       & " (Npix = 12*Nside**2, where Nside HAS to be a power of 2, eg: 32, 512, ...)" )
  nsmax = parse_int(handle, 'nsmax', default=32, descr=description)
  if (nside2npix(nsmax) < 0) then
     print *, " Error: nsmax is not a power of two."
     if (handle%interactive) goto 3
     stop 1
  endif

  !     --- gets the L range for the simulation ---
  WRITE(chline,"(a,i5,a)") "We recommend: (0 <= l <= l_max <= ",3*nsmax-1,")"
  description = concatnl(&
       & "", &
       & " Enter the maximum l range (l_max) for the computation. ", &
       & chline )
  nlmax = parse_int(handle, 'nlmax', default=2*nsmax, descr=description)

  !     --- gets the output sky map filename ---
  description = concatnl(&
       & "", &
       & " Enter Output file name (eg, plm.fits) :", &
       & "  (or !plm.fits to overwrite an existing file)" )
  outfile = parse_string(handle, "outfile", &
       & default="!plm.fits", descr=description, filestatus="new")

  call parse_check_unused(handle, code=lcode)
  call parse_summarize(handle, code=lcode)
  call parse_finish(handle)

  nmmax = nlmax ! not required
  n_plm_8 = (nmmax+1_i8b)*(2*nlmax-nmmax+2_i8b)*nsmax
  n_plm   = (nmmax+1    )*(2*nlmax-nmmax+2    )*nsmax

!   if (n_plm_8 > 2147483647 .or. n_plm < 0) then
!      print*,' '
!      write(*,'(a,i6,i6,a)') 'Error: these values of Nside and l_max (',nsmax,nlmax,')'
!      print*,' are too large for the current implementation of '//trim(lcode)//'.'
!      print*,'See documentation for details.'
!      call fatal_error('Aborting')
!   endif
  if (n_plm_8 > 2147483647 .or. n_plm < 0) then
     print*,' '
     write(*,'(a,i6,i6,a)') 'Warning: these values of Nside and l_max (',nsmax,nlmax,')'
     print*,' are too large for the current implementation of '//trim(lcode)//'.'
     print*,'See documentation for details.'
     print*,'This may create problems.'
  endif

  nd = 1
  if (polarisation) nd = 3
  write(*,'(a,f10.2,a)') 'The code requires ',n_plm_8 * nd * 8e-6, ' Mbytes'

  ALLOCATE(plm(0:n_plm_8-1,1:nd),stat = status)
  call assert_alloc(status,code,'plm')
  ! print*,'allocated'

  call plm_gen(nsmax, nlmax, nmmax, plm)

  header = ' '
  call add_card(header,'CREATOR', code,        'Software creating the FITS file')
  call add_card(header,'VERSION', version,     'Version of the simulation software')
  ! NSIDE, MAX-LPOL, MAX-MPOL, POLAR kwds are automatically updated (or added) by write_plm
  call add_card(header,'NSIDE'   ,nsmax,   'Resolution parameter for HEALPIX')
  call add_card(header,'MAX-LPOL',nlmax      ,'Maximum multipole l of P_lm')
  call add_card(header,'MAX-MPOL',nmmax      ,'Maximum m of P_lm')
  call add_card(header,'POLAR',   polarisation, 'Polarisation included (T/F)')
  call add_card(header)

  call add_card(header,'COMMENT','-----------------------------------------------')
  call add_card(header,'COMMENT','     Data Description Specific Keywords       ')
  call add_card(header,'COMMENT','-----------------------------------------------')
  call add_card(header)
!   call add_card(header,'COMMENT','contains precomputed Y_lms,')
!   call add_card(header,'COMMENT','in the order: ')
!   call add_card(header,'COMMENT','m=0:      Y_00, Y_10, Y_20, Y_30, ...')
!   call add_card(header,'COMMENT','m=1:            Y_11, Y_21, Y_31, ...')
!   call add_card(header,'COMMENT','m=2:                  Y_22, Y_32, ...')
!   call add_card(header)
!   call add_card(header,'COMMENT',' for each ring from North Pole to Equator.')
!   call add_card(header,'COMMENT','synfast and anafast execute faster with precomputed Y_lms')
!   call add_card(header)
!   call add_card(header,'TTYPE1', 'SCALAR',  'Scalar Spherical Harmonics (Y_lm)')
!   if (polarisation) then
!      call add_card(header,'TTYPE2', 'SPINSUM', 'half sum        of Spin-weighted Spherical Harmonics')
!      call add_card(header,'TTYPE3', 'SPINDIFF','half difference of Spin-weighted Spherical Harmonics')
!      call add_card(header,'COMMENT','spin-weighted spherical harmocis now have CMBFAST normalisation')
!   endif
!   call add_card(header)
  call add_card(header,'COMMENT','The complex spherical harmonics are defined as')
  call add_card(header,'COMMENT','  Y_lm(theta, phi) = L_lm(theta) * exp(i m phi),')
  call add_card(header,'COMMENT',' where L_lm are real numbers related to the Legendre Polynomials Plm')
  call add_card(header)
  call add_card(header,'COMMENT','This file contains precomputed L_lms,')
  call add_card(header,'COMMENT','in the order: ')
  call add_card(header,'COMMENT','m=0:      L_00, L_10, L_20, L_30, ...')
  call add_card(header,'COMMENT','m=1:            L_11, L_21, L_31, ...')
  call add_card(header,'COMMENT','m=2:                  L_22, L_32, ...')
  call add_card(header)
  call add_card(header,'COMMENT',' for each ring from North Pole to Equator.')
  call add_card(header,'COMMENT','synfast and anafast execute faster with precomputed L_lms')
  call add_card(header)
  call add_card(header,'TTYPE1', 'SCALAR',  'Scalar Spherical Harmonics (L_lm)')
  if (polarisation) then
     call add_card(header,'TTYPE2', 'SPINSUM', 'half sum   of Spin-weighted Spherical Harmonics')
     call add_card(header,'TTYPE3', 'SPINDIFF','half diff. of Spin-weighted Spherical Harmonics')
     call add_card(header,'COMMENT','spin-weighted spherical harmonics now have CMBFAST normalisation')
  endif
  call add_card(header)


  nlheader = SIZE(header)
  call write_plm(plm, n_plm_8, nd, header, nlheader, outfile, nsmax, nlmax)

  deallocate(plm)

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

  write(*,9005) " Resolution parameter   : ", nsmax
  write(*,9005) " Multipole range        : 0 < l <= ", nlmax
  write(*,9000) " Output file            : "//TRIM(outfile)
  write(*,9030) " Clock and CPU time [s] : ", clock_time, ptime

9000 format(a)
9005 format(a,i8)
9010 format(a,i16)
9020 format(a,g20.5)
9030 format(a,f11.2,f11.2)

end program plmgen
