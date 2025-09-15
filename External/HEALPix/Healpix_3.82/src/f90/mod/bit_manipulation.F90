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
!! This module provides functions that permute bits of the binary
!! representation of a number.
!!
!! @author Benjamin D. Wandelt October 1997
!! edited by E. Hivon, October 2001 to be 'F' compatible
! 2009-06-16: deals with Nside > 8192
module bit_manipulation

use healpix_types
implicit none
private

integer(i4b), parameter :: oddbits=89478485   ! 2^0 + 2^2 + 2^4+..+2^26
integer(i4b), parameter :: evenbits=178956970 ! 2^1 + 2^3 + 2^4+..+2^27

#ifndef NO64BITS
interface swapLSBMSB
   module procedure swapLSBMSB, swapLSBMSB8
end interface
interface invswapLSBMSB
   module procedure invswapLSBMSB, invswapLSBMSB8
end interface
interface invLSB
   module procedure invLSB, invLSB8 
end interface
interface invMSB
   module procedure invMSB, invMSB8
end interface
#endif

public :: swapLSBMSB, invswapLSBMSB, invLSB, invMSB

contains

!! Returns i with even and odd bit positions interchanged.
function swapLSBMSB(i)
  integer(i4b) :: swapLSBMSB
  integer(i4b), intent(in) :: i

  swapLSBMSB = IAND(i,evenbits)/2 + IAND(i,oddbits)*2
end function swapLSBMSB

!! Returns NOT(i) with even and odd bit positions interchanged.
function invswapLSBMSB(i)
  integer(i4b) :: invswapLSBMSB
  integer(i4b), intent(in) :: i

  invswapLSBMSB = NOT(swapLSBMSB(i))
end function invswapLSBMSB

!! Returns i with odd (1,3,5,...) bits inverted.
function invLSB(i)
  integer(i4b) :: invLSB
  integer(i4b), intent(in) :: i

  invLSB = IEOR(i,oddbits)
end function invLSB

!! Returns i with even (0,2,4,...) bits inverted.
function invMSB(i)
  integer(i4b) :: invMSB
  integer(i4b), intent(in) :: i

  invMSB = IEOR(i,evenbits)
end function invMSB

!**********************************************************************
#ifndef NO64BITS
!! Returns i with even and odd bit positions interchanged.
function swapLSBMSB8(i)
  integer(i8b) :: swapLSBMSB8
  integer(i8b), intent(in) :: i
  integer(i8b), parameter :: oddbits8 = 96076792050570581_i8b ! 2^0 + 2^2 + 2^4+..+2^56
  integer(i8b), parameter :: evenbits8=192153584101141162_i8b ! 2^1 + 2^3 + 2^4+..+2^57

  swapLSBMSB8 = IAND(i,evenbits8)/2_i8b + IAND(i,oddbits8)*2_i8b
end function swapLSBMSB8

!! Returns NOT(i) with even and odd bit positions interchanged.
function invswapLSBMSB8(i)
  integer(i8b) :: invswapLSBMSB8
  integer(i8b), intent(in) :: i

  invswapLSBMSB8 = NOT(swapLSBMSB8(i))
end function invswapLSBMSB8

!! Returns i with odd (1,3,5,...) bits inverted.
function invLSB8(i)
  integer(i8b) :: invLSB8
  integer(i8b), intent(in) :: i
  integer(i8b), parameter :: oddbits8 = 96076792050570581_i8b ! 2^0 + 2^2 + 2^4+..+2^56

  invLSB8 = IEOR(i,oddbits8)
end function invLSB8

!! Returns i with even (0,2,4,...) bits inverted.
function invMSB8(i)
  integer(i8b) :: invMSB8
  integer(i8b), intent(in) :: i
  integer(i8b), parameter :: evenbits8=192153584101141162_i8b ! 2^1 + 2^3 + 2^4+..+2^57

  invMSB8 = IEOR(i,evenbits8)
end function invMSB8
#endif
!**********************************************************************

end module bit_manipulation
