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
module synmod

  USE healpix_types
  USE alm_tools, only: create_alm, alm2map, alm2map_der, pixel_window, generate_beam
  USE pix_tools, only: nside2npix
  USE fitstools, only: write_bintab, number_of_alms, read_par, read_conbintab, read_bintab, dump_alms, getsize_fits
  USE head_fits, only: add_card, merge_headers, del_card, write_minimal_header
  use rngmod,    only: rand_init, rand_gauss, planck_rng
  USE misc_utils,only: assert_alloc, fatal_error, wall_clock_time, brag_openmp, string
  USE paramfile_io, only: paramfile_handle, parse_init, parse_int, &
         parse_string, parse_double, parse_lgt, parse_summarize, parse_check_unused, &
         parse_finish, concatnl, scan_directories, get_healpix_data_dir, get_healpix_test_dir, &
         get_healpix_pixel_window_file
  implicit none

  character(len=*), parameter :: CODE = "SYNFAST"
  character(len=FILENAMELEN)  :: lcode

  private
  public :: syn_sub_s, syn_sub_d, CODE, lcode

  contains

    subroutine syn_sub_s(parafile)
      ! single precision implementation
      integer(I4B),               parameter  :: KMAP  = SP  ! precision of map arrays
      character(len=FILENAMELEN), intent(in) :: parafile
      include 'syn_sub_inc.f90'
    end subroutine syn_sub_s
    !------------------------------------


    subroutine syn_sub_d(parafile)
      ! double precision implementation
      integer(I4B),               parameter  :: KMAP  = DP  ! precision of map arrays
      character(len=FILENAMELEN), intent(in) :: parafile
      include 'syn_sub_inc.f90'
    end subroutine syn_sub_d

end module synmod
!====================================================================================

program synfast

  use healpix_types
  use extension, only: getArgument, nArguments
  use synmod,    only: syn_sub_s, syn_sub_d, CODE, lcode
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
     call syn_sub_d(parafile)
  else
     call syn_sub_s(parafile)
  endif

  stop
end program synfast
