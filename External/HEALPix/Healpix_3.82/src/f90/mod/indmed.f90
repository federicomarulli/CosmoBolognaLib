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
!
! module m_indmed of Orderpack 2.0, by Michel Olagnon   http://www.fortran-2000.com/rank/
!
! edited 2012-07-16: 
!         - use includes to avoid repetition of same lines of code 
!           for R,D, I (and forthcoming) implementations
!         - remove SAVE attribute of IDONT, and remove IDONT from module header
!
! 2016-10-05: faster code by turning IDONT array from allocatable to stack array
!
Module m_indmed
  Integer, Parameter :: kdp = selected_real_kind(12,200)
  Integer, Parameter :: ksp = selected_real_kind(5, 30)
  Integer, Parameter :: i4b = selected_int_kind(9)
  public :: indmed
  private :: kdp, ksp, i4b
  private :: R_indmed, I_indmed, D_indmed
  private :: r_med, i_med, d_med
!!  Integer, Allocatable, Dimension(:), Private, Save :: IDONT
  interface indmed
     module procedure d_indmed, r_indmed, i_indmed
  end interface
contains

! *************************************************************
  Subroutine D_indmed (XDONT, INDM)
    !  Returns index of median value of XDONT.
    ! __________________________________________________________
    Real (kind=kdp), Dimension (:), Intent (In) :: XDONT
    Integer (kind=i4b),             Intent (Out):: INDM
    integer, dimension(size(xdont)) :: idont
    Integer :: tmpout, IDON ! can be 4 or 8 bytes depending on default
    ! __________________________________________________________
    !
    Do IDON = 1, SIZE(XDONT)
       IDONT (IDON) = IDON
    End Do
    !
    Call d_med (XDONT, IDONT, tmpout)
    indm = tmpout
    !
  End Subroutine D_indmed

  Recursive Subroutine d_med (XDATT, IDATT, ires_med)
    !  Finds the index of the median of XDONT using the recursive procedure
    !  described in Knuth, The Art of Computer Programming,
    !  vol. 3, 5.3.3 - This procedure is linear in time, and
    !  does not require to be able to interpolate in the
    !  set as the one used in INDNTH. It also has better worst
    !  case behavior than INDNTH, but is about 30% slower in
    !  average for random uniformly distributed values.
    ! __________________________________________________________
    Real (kind=kdp), Dimension (:), Intent (In) :: XDATT
    Integer, Dimension (:), Intent (In) :: IDATT
    Integer, Intent (Out):: ires_med
    ! __________________________________________________________
    !
    Real (kind=kdp), Parameter :: XHUGE = HUGE (XDATT)
    Real (kind=kdp) :: XWRK, XWRK1, XMED7, XMAX, XMIN
    !
    Integer, Dimension (7*(((Size (IDATT)+6)/7+6)/7)) :: ISTRT, IENDT, IMEDT
    Integer, Dimension (7*((Size(IDATT)+6)/7)) :: IWRKT
    Integer :: NTRI, NMED, NORD, NEQU, NLEQ, IMED, IDON, IDON1
    Integer :: IDEB, ITMP, IDCR, ICRS, ICRS1, ICRS2, IMAX, IMIN
    Integer :: IWRK, IWRK1, IMED1, IMED7, NDAT
    !
    NDAT = Size (IDATT)
    include 'indmed_part1.f90'
    Call d_med (XDATT, IMEDT(1:IDON1), IMED7)
    include 'indmed_part2.f90'
    !
  END Subroutine d_med

  ! *************************************************************
  Subroutine R_indmed (XDONT, INDM)
    !  Returns index of median value of XDONT.
    ! __________________________________________________________
    Real (kind=ksp), Dimension (:), Intent (In) :: XDONT
    Integer (kind=i4b),             Intent (Out):: INDM
    integer, dimension(size(xdont)) :: idont
    Integer :: tmpout, IDON ! can be 4 or 8 bytes depending on default
    ! __________________________________________________________
    !
    Do IDON = 1, SIZE(XDONT)
       IDONT (IDON) = IDON
    End Do
    !
    Call r_med (XDONT, IDONT, tmpout)
    indm = tmpout
    !
  End Subroutine R_indmed
  Recursive Subroutine r_med (XDATT, IDATT, ires_med)
    ! __________________________________________________________
    Real, Dimension (:), Intent (In) :: XDATT
    Integer, Dimension (:), Intent (In) :: IDATT
    Integer, Intent (Out) :: ires_med
    ! __________________________________________________________
    !
    Real, Parameter :: XHUGE = HUGE (XDATT)
    Real :: XWRK, XWRK1, XMED7, XMAX, XMIN
    !
    Integer, Dimension (7*(((Size (IDATT)+6)/7+6)/7)) :: ISTRT, IENDT, IMEDT
    Integer, Dimension (7*((Size(IDATT)+6)/7)) :: IWRKT
    Integer :: NTRI, NMED, NORD, NEQU, NLEQ, IMED, IDON, IDON1
    Integer :: IDEB, ITMP, IDCR, ICRS, ICRS1, ICRS2, IMAX, IMIN
    Integer :: IWRK, IWRK1, IMED1, IMED7, NDAT
    !
    NDAT = Size (IDATT)
    include 'indmed_part1.f90'
    Call r_med (XDATT, IMEDT(1:IDON1), IMED7)
    include 'indmed_part2.f90'
    !
  END Subroutine r_med

  ! *************************************************************
  Subroutine I_indmed (XDONT, INDM)
    !  Returns index of median value of XDONT.
    ! __________________________________________________________
    Integer, Dimension (:), Intent (In) :: XDONT
    Integer, Intent (Out) :: INDM
    integer, dimension(size(xdont)) :: idont
    Integer :: tmpout, IDON
    ! __________________________________________________________
    !
    Do IDON = 1, SIZE(XDONT)
       IDONT (IDON) = IDON
    End Do
    !
    Call i_med (XDONT, IDONT, tmpout)
    INDM = tmpout
!
  End Subroutine I_indmed

  Recursive Subroutine i_med (XDATT, IDATT, ires_med)
    ! __________________________________________________________
    Integer, Dimension (:), Intent (In) :: XDATT
    Integer, Dimension (:), Intent (In) :: IDATT
    Integer, Intent (Out) :: ires_med
    ! __________________________________________________________
    !
    Integer, Parameter :: XHUGE = HUGE (XDATT)
    Integer :: XWRK, XWRK1, XMED7, XMAX, XMIN
    !
    Integer, Dimension (7*(((Size (IDATT)+6)/7+6)/7)) :: ISTRT, IENDT, IMEDT
    Integer, Dimension (7*((Size(IDATT)+6)/7)) :: IWRKT
    Integer :: NTRI, NMED, NORD, NEQU, NLEQ, IMED, IDON, IDON1
    Integer :: IDEB, ITMP, IDCR, ICRS, ICRS1, ICRS2, IMAX, IMIN
    Integer :: IWRK, IWRK1, IMED1, IMED7, NDAT
    !
    NDAT = Size (IDATT)
    include 'indmed_part1.f90'
    Call i_med (XDATT, IMEDT(1:IDON1), IMED7)
    include 'indmed_part2.f90'
  END Subroutine i_med

end module m_indmed
