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
module anamod

  USE healpix_types
  USE alm_tools,  only: map2alm_iterative, alm2cl
  USE fitstools,  only: getsize_fits, input_map, read_par, read_bintab, write_asctab, dump_alms, &
       check_input_map, unfold_weightsfile !, unfold_weightslist
  USE head_fits,  only: add_card, write_minimal_header
  use misc_utils, only: assert_alloc, fatal_error, wall_clock_time, brag_openmp, string
  USE pix_tools,  only: convert_nest2ring, convert_ring2nest, nside2npix, remove_dipole, vec2ang
  USE paramfile_io, only: paramfile_handle, parse_init, parse_int, &
         parse_string, parse_double, parse_summarize, parse_check_unused, &
         parse_finish, concatnl, scan_directories, get_healpix_data_dir, get_healpix_test_dir, &
         get_healpix_weight_file
  use udgrade_nr, only: udgrade_nest, udgrade_ring
  implicit none

  character(len=*), parameter :: CODE = "ANAFAST"
  character(len=FILENAMELEN)  :: lcode

  private
  public :: ana_sub_s, ana_sub_d, CODE, lcode, check_input_map

  contains

    subroutine ana_sub_s(parafile)
      ! single precision implementation
      integer(I4B),               parameter  :: KMAP  = SP  ! precision of map arrays
      character(len=FILENAMELEN), intent(in) :: parafile
      include 'ana_sub_inc.f90'
    end subroutine ana_sub_s
    !------------------------------------


    subroutine ana_sub_d(parafile)
      ! double precision implementation
      integer(I4B),               parameter  :: KMAP  = DP  ! precision of map arrays
      character(len=FILENAMELEN), intent(in) :: parafile
      include 'ana_sub_inc.f90'
    end subroutine ana_sub_d
    !------------------------------------



end module anamod
!====================================================================================

program anafast

  use healpix_types
  use extension, only: getArgument, nArguments
  use anamod,    only: ana_sub_s, ana_sub_d, CODE, lcode, check_input_map
  use misc_utils,only: assert, fatal_error, strlowcase

  integer(i4b)                :: n_args, i
  character(len=FILENAMELEN)  :: arg
  logical(lgt)                :: do_double
  character(len=FILENAMELEN)  :: parafile

  ! count arguments, should be 0, 1 or 2
  n_args = nArguments()
  lcode = trim(strlowcase(CODE))
  call assert(n_args <= 2,' Usage: '//trim(lcode)//' [-s|--single|-d|--double] [parameter_file_name]')

  parafile = ''
  do_double = .false.

  ! parse arguments
  do i=1, n_args
     call getArgument(i, arg)
     arg = trim(adjustl(arg))
     if (arg(1:1) == '-') then
        select case (trim(arg))
        case ('-d', '--double')
           do_double = .true.
        case ('-s', '--single')
           do_double = .false.
        case default
           call fatal_error(CODE//': Invalid argument: '//trim(arg))
        end select
     else
        parafile = arg
     endif
  enddo

  ! start calculations
  if (do_double) then
     call ana_sub_d(parafile)
  else            
     call ana_sub_s(parafile)
  endif

  stop
end program anafast
