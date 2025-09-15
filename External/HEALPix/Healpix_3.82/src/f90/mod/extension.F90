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
module extension
  !
  ! defines in F90 some C commands
  ! These extensions are not completely standard in all F90/95 compilers
  !
  ! getEnvironment   : emulates getenv
  ! getArgument      : emulates getarg
  ! nArguments       : emulates iargc
  !
  ! written by Eric Hivon, Nov 2001
  !
  ! exit_with_status : verbose and clean exit, added by M.R.
  ! 2005-08: edited for Gfortran
  ! 2013-05-07: G95-compatible
  ! 2015-07-31: G95-compatible
  ! 2016-05: edited for __GFORTRAN__
  ! 2019-02: M.R. replaced C with Fortran2003 commands

  USE healpix_types, ONLY : I4B, I8B
  IMPLICIT none

  private
  public :: getEnvironment, getArgument, nArguments, exit_with_status

  ! interface for the exit() method, which is part of the C library and
  ! should be linked automatically to Fortran programs
  interface
    subroutine exit(status) bind(C, name="exit")
    use iso_c_binding
    integer(c_int), value :: status
    end subroutine
  end interface

  contains

    ! ===========================================================
    function nArguments() result(narg)
      ! ===========================================================
      integer(kind=I4B) :: narg
      ! ===========================================================
      narg = command_argument_count()
    end function nArguments
    ! ===========================================================
    subroutine getEnvironment(name, value)
      ! ===========================================================
      character(len=*), intent(in) :: name
      character(len=*), intent(out) :: value
      value=""
      CALL get_environment_variable(name, value)
    end subroutine getEnvironment
    ! ===========================================================
    subroutine getArgument(index, argument)
      ! ===========================================================
      integer(kind=I4B), intent(in) :: index
      character(len=*), intent(out) :: argument
      ! ===========================================================
      argument=""
      call get_command_argument(index, argument)
    end subroutine getArgument

    ! ===========================================================
    subroutine exit_with_status (code, msg)
      ! ===========================================================
      integer(i4b), intent(in) :: code
      character (len=*), intent(in), optional :: msg
      ! ===========================================================
      if (present(msg)) print *,trim(msg)
      print *,'program exits with exit code ', code
      call exit (code)
    end subroutine exit_with_status

end module extension
