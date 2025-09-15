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
module altmod

  use healpix_types
  use alm_tools, only: alter_alm, pixel_window, generate_beam, rotate_alm
  use fitstools, only: read_par, number_of_alms, getsize_fits, alms2fits, fits2alms, dump_alms
  use head_fits, only: add_card !, write_minimal_header
  use misc_utils,only: assert, assert_alloc, fatal_error, wall_clock_time, string, strupcase, brag_openmp
  use paramfile_io, only: paramfile_handle, parse_init, parse_int, &
         parse_string, parse_double, parse_summarize, parse_check_unused, parse_finish, concatnl, &
         get_healpix_pixel_window_file, get_healpix_data_dir, scan_directories
  use coord_v_convert, only: coordsys2euler_zyz
  implicit none

  character(len=*), parameter :: CODE = "ALTERALM"
  character(len=FILENAMELEN)  :: lcode

!   interface alt_sub
!      module procedure alt_sub_s, alt_sub_d
!   end interface

  private
  public alt_sub_s, alt_sub_d, CODE, lcode

  contains

    subroutine alt_sub_s(parafile)
      ! single precision implementation
      integer(I4B),               parameter  :: KMAP  = SP  ! precision of map arrays
      character(len=FILENAMELEN), intent(in) :: parafile
      include 'alt_sub_inc.f90'
    end subroutine alt_sub_s
    !------------------------------------


    subroutine alt_sub_d(parafile)
      ! double precision implementation
      integer(I4B),               parameter  :: KMAP  = DP  ! precision of map arrays
      character(len=FILENAMELEN), intent(in) :: parafile
      include 'alt_sub_inc.f90'
    end subroutine alt_sub_d

end module altmod
!====================================================================================

program alteralm

  use healpix_types
  use extension, only: getArgument, nArguments
  use altmod,    only: alt_sub_s, alt_sub_d, CODE, lcode
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
        case ('-d')
           do_double = .true.
        case ('--double')
           do_double = .true.
        case ('-s')
           do_double = .false.
        case ('--single')
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
     call alt_sub_d(parafile)
  else
     call alt_sub_s(parafile)
  endif

  stop
end program alteralm
