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
Module sky_sub

  !Subroutines used by sky_ng_sim and sky_ng_sim_bin

Contains
  
  !--------------------------------------------------------
  
  Subroutine add_subscript(text, number)

    Implicit none

    Character(len = 20), Intent(InOut) :: text
    Integer, Intent(In) :: number

    Integer i
    Character(len = 10) :: form

    i = Ceiling(log10(number + 0.5))
    If (i .ge. 10) Stop "A ridiculously large number of bins!"

    Write (form, "(a,i1,a)") "(a,a,i",i,")"

    Write (text, form) TRIM(text), "_", number

  End Subroutine add_subscript

  !--------------------------------------------------------

  Subroutine read_powerspec(infile, nsmax, nlmax, cls, header_PS, fwhm, units_power, winfile, beam_file)

    Use healpix_types
    USE fitstools, ONLY : fits2cl
    Use alm_tools, Only : generate_beam, pixel_window
    use utilities, only : die_alloc
    Use head_fits, only : add_card, get_card, merge_headers

    Implicit none

    Character(len = *), Intent(In) :: infile
    Integer(I4B), Intent(In)       :: nsmax
    Integer(I4B), Intent(InOut)    :: nlmax
    Real(SP), dimension(0:nlmax,1:1), Intent(Out) :: cls
    Character(Len=80), dimension(1:), Intent(Out) :: header_PS
    Real(DP) :: fwhm
    Character(len=80), dimension(1:1), Intent(Out) :: units_power
    Character(len=*), Intent(In), optional :: winfile
    Character(len=*), Intent(In), optional :: beam_file


    Logical fitscl
    Integer(I4B) :: nlheader
    Integer(I4B) :: l_max, npw
    Real(DP) :: lbeam
    Integer i 
    Integer(I4B) :: status, count_tt = 0
    Real(SP) :: quadrupole
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: pixlw, beamw
    CHARACTER(LEN=*), PARAMETER :: code = 'READ_POWERSPEC'
    CHARACTER(LEN=80), DIMENSION(1:180)          :: header_file
    CHARACTER(LEN=20)                            ::  string_quad
    !  CHARACTER(LEN=80), DIMENSION(:), allocatable :: units_power
    CHARACTER(LEN=80) :: temptype, com_tt

    ! set maximum multipole (depends on user choice and pixel window function)
    ! if using pixel window function l_max can't be greater than 4*nsmax
    l_max = nlmax
    npw = 4*nsmax + 1
    if (present(winfile)) then
       l_max = min(l_max, npw - 1)
    else
       npw = l_max + 1
    endif
    if (l_max < nlmax) then
       print*,'C_l are only non-zero for 0 < l <= ',l_max
    endif
    nlmax = l_max

    !-----------------------------------------------------------------------

    cls = 0.0
    nlheader = SIZE(header_file)
    header_file = ''
    ! Find out whether input cl file is a fits file
    fitscl = (index(infile, '.fits') /= 0)
    If (fitscl) Then
!!!       call read_asctab(infile, cls, nlmax, 1, header_file, nlheader, units=units_power)
       call fits2cl(infile, cls, nlmax, 1, header_file, units=units_power)
       call get_card(header_file,"TEMPTYPE",temptype,com_tt,count=count_tt)
    Else
       call read_camb(infile, cls, nlmax)
       units_power = 'microKelvin^2'
       !Check to see if enough values were present in file
       If (nlmax .LT. l_max) then
          Write (*,*) "Maximum l-value in cl file is ",nlmax
          Write (*,*) "Power spectrum above this value set to zero"
          l_max = nlmax
       end if
    End If

    if (count_tt < 1) then
       temptype = "THERMO"
       com_tt = " temperature : either THERMO or ANTENNA"
       print*,'   Will assume '//temptype
    endif

    quadrupole = cls(2,1)
    !     --- creates the header relative to power spectrum ---
    header_PS = ''
    call add_card(header_PS)
    call add_card(header_PS,'HISTORY',' alm generated from following power spectrum')
    call add_card(header_PS)
    call add_card(header_PS,'COMMENT','----------------------------------------------------')
    call add_card(header_PS,'COMMENT','Planck Power Spectrum Description Specific Keywords')
    call add_card(header_PS,'COMMENT','----------------------------------------------------')
    call add_card(header_PS,'COMMENT','Input power spectrum in : ')
    call add_card(header_PS,'COMMENT',TRIM(infile))
    call add_card(header_PS)
    call add_card(header_PS,'COMMENT','Quadrupole')
    write(string_quad,'(1pe15.6)') quadrupole
    call add_card(header_PS,'COMMENT','  C(2) = '//string_quad//' '//trim(units_power(1)))
    call add_card(header_PS)
    call add_card(header_PS,"TEMPTYPE",temptype,com_tt)
    call add_card(header_PS)
    ! ------- insert header read in power spectrum file -------
    call merge_headers(header_file, header_PS)
    !    deallocate(units_power)

    ! define beam
    allocate(beamw(0:l_max,1:1), stat = status)
    if (status /= 0) call die_alloc(code,'beamw')
    If ((fwhm .NE. 0d0) .OR. present(beam_file)) Then
       call generate_beam(real(fwhm,kind=dp), l_max, beamw, beam_file)
    Else
       beamw(:,:) = 1.0_dp
    End If

    ! get the pixel window function 
    allocate(pixlw(0:npw-1,1:1), stat = status)
    if (status /= 0) call die_alloc(code,'pixlw')
    if(present(winfile)) then
       call pixel_window(pixlw, windowfile=winfile)
    else
       pixlw(:,:)  = 1.0_dp
    endif

    ! Multiply by beam and pixel window functions
!    cls(0:l_max, 1) = cls(0:l_max,1)*beamw(0:l_max,1)*pixlw(0:l_max,1) !EH-2008-03-05
    cls(0:l_max, 1) = cls(0:l_max,1)* (beamw(0:l_max,1)*pixlw(0:l_max,1))**2

    deallocate(beamw)
    deallocate(pixlw)

  End Subroutine read_powerspec

  !--------------------------------------------------------

  Subroutine read_camb(infile, cl_array, lmax)

    Use healpix_types

    Implicit none

    Character (len = *), intent(IN) :: infile
    Integer(I4B), intent(INOUT) :: lmax
    Real(SP), Intent(OUT), Dimension(0:lmax,1:1) :: cl_array

    Real(DP) :: clval
    Integer :: lval = 0
    Real(DP) :: cmbT = 2.726d0

    ! Read in the cl values and convert to units of uK^2
    Open (Unit = 23, file = infile, status = 'old', action = 'read', err = 9980)
    Do while (lval .LT. lmax)
       Read (23, *, end = 800) lval, clval
       cl_array(lval, 1) = TWOPI*clval/(lval*(lval+1))*cmbT*cmbT*1.d12
    End Do
800 lmax = lval

    Close (23)
    Return

9980 Stop 'Error opening cl file to read'

  End Subroutine read_camb


End Module sky_sub
