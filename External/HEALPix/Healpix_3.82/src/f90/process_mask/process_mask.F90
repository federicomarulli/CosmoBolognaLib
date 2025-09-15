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
!=======================================================================
!
! For a given input binary mask, in which pixels have either value 0 (=invalid) or 1 (=valid),
! this code produces a map containing for each valid pixel,
! its distance (in Radians) to the closest invalid pixel.
! It can then be used to define an apodisation scheme.
! Distances are measured between pixel centers.
!
! It is possible to treat small holes (=cluster of adjacent invalid pixels) as valid,
! by specifying a minimal number of pixels and/or minimal surface area (which ever is the largest),
! and the resulting new mask can be output.
!
! The output FITS files have the same ordering as the input mask
! (even though the processing is done in NESTED ordering).
!
! 2 pixels are adjacent if they have at least one point in common 
! (eg, a pixel corner or a pixel side).
!
! 
! HISTORY:
!  E. Hivon 2006-2010 firt developments
!  2011-11-22: inclusion in Healpix 
!=======================================================================
program process_mask

  use healpix_types
  use extension, only: getArgument, nArguments
  use fitstools, only: getsize_fits, input_map, output_map
  use head_fits, only: add_card, write_minimal_header
  use pix_tools, only: convert_ring2nest, convert_nest2ring, nside2npix
  USE paramfile_io, only: paramfile_handle, parse_init, parse_int, &
         parse_string, parse_double, parse_summarize, parse_check_unused, parse_finish, concatnl
  use misc_utils, only: assert_alloc, fatal_error, wall_clock_time, brag_openmp
  use mask_tools, only: fill_holes_nest, dist2holes_nest
  implicit none

  integer(I4B) :: nbordpix, status, n_args
  integer(I4B) :: nside, npix, order, new_min_size
  integer(I4B), allocatable, dimension(:) :: mask
  real(SP), allocatable, dimension(:,:) :: map
  real(DP), allocatable, dimension(:,:) :: map2
  real(DP):: surface_deg2, surf_pix_deg
  
  real(kind=SP)                 :: clock_time, time0, time1, t0, t1
  real(kind=SP)                 :: ptime, ptime0, ptime1, pt0, pt1
  type(paramfile_handle)     :: handle

  character(len=FILENAMELEN) :: parafile
  character(len=FILENAMELEN)          :: description
  character(len=80), dimension(1:200) :: header
  character(len=*), parameter :: code = 'process_mask'
  character(len=*), parameter :: VERSION = HEALPIX_VERSION
  integer(i8b) :: tot1, tot2
  
  type proc_mask
     character(len=filenamelen) :: mask_file, distance_file, filled_file !, file_bord
     real(DP)                   :: hole_min_surf_arcmin2 
     integer(i4b)               :: hole_min_size
  end type proc_mask
  type(proc_mask) :: pm

!------------------------------------------------
  call wall_clock_time(time0)
  call cpu_time(ptime0)

  print *, " "
  print *, "                "//code//" "//version
  write(*,'(a)') " *** Mask processing ***"
  ! ------ parse arguments -----
  n_args = nArguments()
  parafile = ''
  if (n_args > 1) then
     stop 'Expected 0 or 1 argument'
  endif
  if (n_args == 1) then
     call getArgument(1, parafile)
  endif
  handle = parse_init(parafile)
  
  description = concatnl('',&
       & 'Enter the input (binary) mask file')
  pm%mask_file = parse_string(handle, 'mask_file', descr=description, filestatus='old')

  description = concatnl('',&
       & 'Enter the minimal size (in pixels) of invalid regions to be kept',&
       & '(can be used together with hole_min_surf_arcmin2 below,'// &
       & ' the result will be the largest of the two)')
  pm%hole_min_size = parse_int(handle,'hole_min_size',descr=description, default=0)

  description = concatnl('',&
       & 'Enter the minimal surface area (in arcmin^2) of invalid regions to be kept',&
       & '(can be used together with hole_min_size above,' // &
       & ' the result will be the largest of the two)')
  !surface_deg2 = parse_double(handle,'surface_deg2',descr=description, default=0.d0)
  pm%hole_min_surf_arcmin2 = parse_double(handle,'hole_min_surf_arcmin2',descr=description, default=0.d0)

  description = concatnl('',&
       & 'Enter the OPTIONAL output file to contain mask with filled-in small invalid regions',&
       & '(as defined above)')
  pm%filled_file = parse_string(handle,'filled_file',descr=description,filestatus='new',default="''")

  description = concatnl('',&
       & 'Enter the OPTIONAL output distance file')
  pm%distance_file = parse_string(handle,'distance_file',descr=description,filestatus='new',default="''")

  !call parse_check_unused(handle, code=code)
  call parse_summarize(handle,code=code)
  call parse_finish(handle)
  ! ----------------------------
  call brag_openmp()
  ! input mask
  npix = getsize_fits(pm%mask_file, nside=nside, ordering=order)
  allocate(map(0:npix-1,1:1),stat=status)
  call assert_alloc(status,code,'map')
  call input_map(pm%mask_file, map, npix, 1)

  surf_pix_deg = fourpi/npix * rad2deg**2 ! area/pixel [Deg^2]
  surface_deg2 = pm%hole_min_surf_arcmin2 / 3600.0_dp
  new_min_size = max( NINT(surface_deg2/surf_pix_deg), pm%hole_min_size)
  ! print*,'hole minimal size (in pixels) :',new_min_size

  if (order == 1) then
     write(*,*) code//'> Converting RING -> NESTED'
     call convert_ring2nest(nside, map)
     !print*,'conversion done'
  endif
  allocate(mask(0:npix-1),stat=status)
  call assert_alloc(status,code,'mask')
  mask = NINT(map(:,1))
  tot1 = sum(mask)

  !----------------------------------------------------------------------
  ! generate new mask by filling small holes
  !----------------------------------------------------------------------
  call fill_holes_nest(nside, new_min_size, mask, mask)
  tot2 = sum(mask)

  if (trim(pm%filled_file) /= '') then
     ! output new mask
     map(:,1) = mask
     if (order == 1) call convert_nest2ring(nside, map)
     
     !print*, trim(file_filled)
     !header(:) = ' '
     call write_minimal_header(header,'MAP',&
          &   nside=nside, &
          &   order=order, &
          &   creator=code, &
          &   version=version)
     call add_card(header,'HISTORY','Hole filler')
     call add_card(header,'HISTORY','Input mask:')
     call add_card(header,'HISTORY',trim(pm%mask_file))
     !call add_card(header,'SURFMIN',surface_deg2,'[Deg^2] Minimal surface')
     call add_card(header,'NP_MIN',new_min_size,'[pixels] Minimal size')
     call output_map(map, header, pm%filled_file)
  endif

  !----------------------------------------------------------------------
  ! compute distance to nearest hole on new mask
  !----------------------------------------------------------------------
  if (trim(pm%distance_file) /= '') then

     allocate(map2(0:npix-1,1:1),stat=status)
     call assert_alloc(status,code,'map2')
     call dist2holes_nest(nside, mask, map2(:,1))

     if (order == 1) call convert_nest2ring(nside, map2)

     call write_minimal_header(header,'MAP',&
          &   nside=nside, &
          &   order=order, &
          &   creator=code, &
          &   version=version, &
          &   units='RADIANS')
     call add_card(header,'TTYPE1','ANGULAR_DISTANCE','[RAD]')
     call add_card(header,'HISTORY','Distance to closest hole for')
     call add_card(header,'HISTORY','Input mask:')
     call add_card(header,'HISTORY',trim(pm%mask_file))
     call add_card(header,'HISTORY','Filled mask:')
     call add_card(header,'HISTORY',trim(pm%filled_file))
     !call add_card(header,'SURFMIN',surface_deg2,'[Deg^2] Minimal surface')
     call add_card(header,'NP_MIN',new_min_size,'[pixels] Minimal size')
     call output_map(map2, header, pm%distance_file)
  endif

  call wall_clock_time(time1)
  call cpu_time(ptime1)
  !----------------------------------------------------------------------
  ! report card
  !----------------------------------------------------------------------
  write(*,9000) " "
  write(*,9000) " Report Card for "//code//" run"
  write(*,9000) "--------------------------------------"
  write(*,9000) " "
  write(*,9000) " Input binary mask       : "//TRIM(pm%mask_file)
  write(*,9010) " Total number of pixels  : ", npix
  write(*,9010) " Input valid pixels      : ", tot1
  write(*,9010) " Effective min. hole size: ", new_min_size
  write(*,9010) " Filled-in pixels        : ", tot2-tot1
  if (trim(pm%filled_file) /= '') write(*,9000) &
       &        " Output filled-in mask    : "//TRIM(pm%filled_file)
  if (trim(pm%distance_file) /= '') write(*,9000) &
       &        " Output distance file    : "//TRIM(pm%distance_file)
  write(*,9030) " Clock and CPU time [s]  : ", time1-time0, ptime1-ptime0

  !-----------------------------------------------------------------------
  !                       end of routine
  !-----------------------------------------------------------------------

  write(*,9000) " "
  write(*,9000) " "//code//"> normal completion"

9000 format(a)
9010 format(a,i16)
9020 format(a,g20.5)
9030 format(a,f11.2,f11.2)



stop
end program process_mask

