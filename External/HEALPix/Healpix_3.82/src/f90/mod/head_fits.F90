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
module head_fits

  !-------------------------------------------------------------------------
  ! this module (head_fits) introduced in Healpix 1.2
  ! merges the module wrap_fits present in version 1.1
  ! and the modules strip_fits that appeared in some intermediate versions
  !
  ! It is designed to read and write header information in FITS files
  !
  ! EH, 2002-08-09
  ! 2007-09-20 added write_minimal_header
  ! 2009-01-08 corrected typo in FITS header written by write_minimal_header (DY_Y -> DY_T)
  ! 2013-05-07: G95-compatible
  ! 2014-05-02: fixed problem with keywords having a long string value
  ! 2018-05-28: write_minimal_header: in agreement with IDL, healpy, C++?
  !    IQU maps are now column names TEMPERATURE, Q_POLARISATION, U_POLARISATION,
  !            
  !-------------------------------------------------------------------------

  ! subroutine add_card   [interface]
  ! subroutine write_hl
  ! subroutine merge_headers
  !
  ! subroutine get_card   [interface]
  !
  ! subroutine del_card
  !
  ! write_minimal_header

  USE healpix_types
  USE misc_utils
  interface add_card
     ! add_card generic routine
     ! ------------------------
     ! adds a card to a fits header
     ! syntax = call add_card(header, keyword, value [, comment])
     !   header  = character string array
     !   keyword = character string (NON BLANK)
     !   value   = can be either logical, integer, real, double, character
     !   comment = (optional) character string scalar
     !
#ifdef NO64BITS
     ! if 64 bit (=8 bytes) integer is NOT supported
     module procedure d_add_card, f_add_card, a_add_card, i_add_card,              l_add_card, v_add_card
#else
     ! if 64 bit (=8 bytes) integer is supported
     module procedure d_add_card, f_add_card, a_add_card, i_add_card, ii_add_card, l_add_card, v_add_card
#endif

  end interface

  interface get_card
     ! get_card generic routine
     ! ------------------------
     ! read a keyword (case UNsensitive) from a FITS header
     ! syntax = call get_card(header, keyword, value [, comment, count])
     !   header  = (input) character string array
     !   keyword = (input) character string scalar
     !   value   = (output) can be either logical, integer, real, double, character
     !              contains on output value of keyword
     !   comment = (output, optional) character string scalar
     !              contains on output value of comment field
     !   count   = (output, optional) integer scalar,
     !       is set to 1 if the keyword is found and 0 otherwise
     !
     ! except if value is declared as a character string,
     ! its type should match that in the fits file
     !
#ifdef NO64BITS
     ! if 64 bit (=8 bytes) integer is NOT supported
     module procedure d_get_card, f_get_card, a_get_card, i_get_card,              l_get_card
#else
     ! if 64 bit (=8 bytes) integer is supported
     module procedure d_get_card, f_get_card, a_get_card, i_get_card, ii_get_card, l_get_card
#endif
  end interface

  interface del_card
     ! del_card generic routine
     ! ------------------------
     ! deletes a card from a fits header
     ! syntax = call del_card(header, keywords) or
     !          call del_card(header, keyword)
     !   header   = character string array
     !   keywords = character string vector (NON BLANK)
     !   keywords = character string  (NON BLANK)
     !
     module procedure del_card1, del_cardn
  end interface

  private

  character (LEN=1 ), private :: dtype
  character (LEN=20), private :: stval
  character (LEN=80), private :: card, stcom
! work arond G95 bug that does not accept SAVEd variables with value 0 or .false. (2013-05-07)
!   integer(kind=I4B) , private :: nlh, statusfits, count_in = 0
!   logical(kind=LGT) , private :: match, exact, casesen = .false., verbose = .true.
  integer(kind=I4B), private :: nlh, statusfits, count_in
  logical(kind=LGT), private :: match, exact
  logical(kind=LGT), private, parameter :: casesen = .false., verbose = .true.

  public :: add_card, merge_headers

  public :: get_card

  public :: del_card

  public :: write_minimal_header
contains

  !=====================================================================
  subroutine d_add_card(header, kwd, value, comment, update) ! double precision
    real(kind=DP), intent(IN) :: value
    character (LEN=80), dimension(:), intent(INOUT) :: header
    character (LEN=*), intent(IN) :: kwd
    character (LEN=*), intent(IN), OPTIONAL :: comment
    logical(LGT),      intent(IN), OPTIONAL  :: update

    write(stval,'(1pe20.12)') value
    call write_hl(header, kwd, stval, comment, update)
    RETURN
  END subroutine d_add_card
  !=====================================================================
  subroutine f_add_card(header, kwd, value, comment, update) ! single precision
    real(kind=SP), intent(IN) :: value
    character (LEN=80), dimension(:), intent(INOUT) :: header
    character (LEN=*), intent(IN) :: kwd
    character (LEN=*), intent(IN), OPTIONAL :: comment
    logical(LGT),      intent(IN), OPTIONAL  :: update

    write(stval,'(1pe20.8)') value
    call write_hl(header, kwd, stval, comment, update)
    RETURN
  END subroutine f_add_card
  !=====================================================================
  subroutine a_add_card(header, kwd, value, comment, update) ! character
    character (LEN=*), intent(IN), OPTIONAL :: value
    character (LEN=80), dimension(:), intent(INOUT) :: header
    character (LEN=*), intent(IN) :: kwd
    character (LEN=*), intent(IN), OPTIONAL :: comment
    logical(LGT),      intent(IN), OPTIONAL  :: update
    character (LEN=240) :: st_value, st_comment
    logical(LGT) :: long_strn, actual

    st_value   = ''
    st_comment = ''
    long_strn=.false.
    actual = (trim(kwd) /= "COMMENT" .and. trim(kwd) /= "HISTORY" .and. trim(kwd) /= "CONTINUE")
    if (present(value)  ) then
       write(st_value,'(a)')   value
       st_value = adjustl(st_value)
       long_strn = (actual .and. len_trim(st_value) > 68)
       if (st_value(1:1) /= "'" .and. actual) st_value = "'"//trim(st_value)//"'"
    endif
    if (present(comment)) write(st_comment,'(a)') comment
    call write_hl(header, kwd, st_value, st_comment, update, long_strn=long_strn)
    RETURN
  END subroutine a_add_card
  !=====================================================================
  subroutine i_add_card(header, kwd, value, comment, update) ! integer (i*4)
    integer(kind=I4B), intent(IN) :: value
    character (LEN=80), dimension(:), intent(INOUT) :: header
    character (LEN=*), intent(IN) :: kwd
    character (LEN=*), intent(IN), OPTIONAL :: comment
    logical(LGT),      intent(IN), OPTIONAL  :: update

    write(stval,'(i20)') value
    call write_hl(header, kwd, stval, comment, update)
    RETURN
  END subroutine i_add_card
  !=====================================================================
#ifndef NO64BITS
  ! compile 8 byte integer flavor only if possible
  subroutine ii_add_card(header, kwd, value, comment, update) ! integer (i*8)
    integer(kind=I8B), intent(IN) :: value
    character (LEN=80), dimension(:), intent(INOUT) :: header
    character (LEN=*), intent(IN) :: kwd
    character (LEN=*), intent(IN), OPTIONAL :: comment
    logical(LGT),      intent(IN), OPTIONAL  :: update

    write(stval,'(i20)') value
    call write_hl(header, kwd, stval, comment, update)
    RETURN
  END subroutine ii_add_card
#endif
  !=====================================================================
  subroutine l_add_card(header, kwd, value, comment, update) ! logical
    logical(kind=LGT), intent(IN) :: value
    character (LEN=80), dimension(:), intent(INOUT) :: header
    character (LEN=*), intent(IN) :: kwd
    character (LEN=*), intent(IN), OPTIONAL :: comment
    logical(LGT),      intent(IN), OPTIONAL  :: update

    write(stval,'(l7)') value
    call write_hl(header, kwd, stval, comment, update)
    RETURN
  END subroutine l_add_card
  !=====================================================================
  subroutine v_add_card(header) ! blank line
    character (LEN=80), dimension(:), intent(INOUT) :: header

    call write_hl(header, 'COMMENT', ' ', ' ')
!    call write_hl(header, ' ', ' ', ' ')
  END subroutine v_add_card
  !=====================================================================
  subroutine write_hl(header, kwd, st_value, comment, update, long_strn)
    IMPLICIT none
    character (LEN=80), dimension(1:), intent(INOUT) :: header
    character (LEN=*), intent(IN) :: kwd
    character (LEN=*), intent(IN), OPTIONAL  :: st_value
    character (LEN=*), intent(IN), OPTIONAL  :: comment
    logical(LGT),      intent(IN), OPTIONAL  :: update, long_strn
    integer(kind=I4B) :: iblank, iwrite, hdtype, kmatch, iwlong

    character (LEN=480) :: fullcard
    character (LEN=80)  :: oldline
    character (LEN=20)  :: kwd2, this_kwd
    logical(LGT)        :: do_update, do_long_strn
    !=====================================================================

    fullcard=''
    oldline=''

    do_update = .false.
    if (present(update)) do_update = update

    do_long_strn = .false.
    if (present(long_strn)) do_long_strn = long_strn

    nlh = size(header)
    ! iblank = first blank line
    iblank = nlh + 1
    do
       if (iblank == 1) exit
       if (trim(header(iblank-1)) /= '') exit
       iblank = iblank - 1
    enddo


    kwd2 = adjustl(strupcase(kwd))
    fullcard = kwd2
    if (present(st_value)) fullcard = fullcard(1:max( 8,len_trim(fullcard)))//' '//trim(adjustl(st_value))
    if (present(comment)) then
       if (len_trim(adjustl(comment)) > 0) then
          if (do_long_strn) then ! put a slash between value and comment for easier idenfication later on
             fullcard = fullcard(1:max(32,len_trim(fullcard)))//' / '//trim(adjustl(comment))
          else
             fullcard = fullcard(1:max(32,len_trim(fullcard)))//' '//trim(adjustl(comment))
          endif
       endif
    endif

    if (iblank > nlh) then
       print*,'WARNING: Header is too short ('//trim(string(nlh))//' lines); card not written:'
       print*,trim(fullcard)
       return
    endif

    if (trim(fullcard) == 'COMMENT') then ! COMMENT alone
       header(iblank) = 'COMMENT'
!       iblank = iblank + 1
       return
    endif

    hdtype = 0
    statusfits = 0
    CALL ftgthd(fullcard(1:79), oldline, hdtype, statusfits)
    if (statusfits == 205) statusfits = 0 ! ignore unmatched quotes due to long strings
    ! hdtype: -1: delete, 0: append or replace, 1: append, 2: END
    select case (hdtype)
    case (-1)! delete keyword
       header(iblank) = fullcard(1:79)
    case (0)
       if (do_update) then ! update (possibly) existing keyword, removing other matches
          kmatch = 0
          do iwrite = 1, iblank ! scan the header
             oldline = header(iwrite)
             this_kwd = strupcase(oldline(1:len_trim(kwd2)))
             if (trim(this_kwd) == trim(kwd2) .or. iwrite == iblank) then ! find first matching line or reach first blank line
                kmatch = kmatch + 1
                if (kmatch == 1) then ! first match: insert or add line
                   call insert_line_in_header(header, iwrite, fullcard, overwrite=do_update, long_strn = do_long_strn)
                else ! other matches : remove header line
                   header(iwrite) = ' '
                   iwlong = iwrite + 1  ! take care of long-string keywords
                   do
                      oldline = adjustl(header(iwlong))
                      if (strupcase(oldline(1:8)) == 'CONTINUE' ) then
                         header(iwlong) = ' '
                         iwlong = iwlong + 1
                      else
                         exit
                      endif
                   enddo

                endif
             endif
          enddo
       else ! append, removing previous matches
          do iwrite = 1, iblank-1
             oldline = header(iwrite)
             this_kwd = strupcase(oldline(1:len_trim(kwd2)))
             if (trim(this_kwd) == trim(kwd2)) then
                header(iwrite) = ' ' ! warning, this does not take care of long-string keywords
                iwlong = iwrite + 1  ! take care of long-string keywords
                do
                   oldline = adjustl(header(iwlong))
                   if (strupcase(oldline(1:8)) == 'CONTINUE' ) then
                      header(iwlong) = ' '
                      iwlong = iwlong + 1
                   else
                      exit
                   endif
                enddo
             endif
          enddo
          call insert_line_in_header(header, iblank, fullcard, overwrite=.false., long_strn = do_long_strn)
       endif
    case (1:2) ! append, no question asked (COMMENT, HISTORY, END)
       call insert_line_in_header(header, iblank, fullcard, overwrite=.false., long_strn = do_long_strn)
    end select
          
    RETURN
  END subroutine write_hl


  !=====================================================================
  subroutine insert_line_in_header(header, idx, card, overwrite, long_strn)
    IMPLICIT none
    character (LEN=80), dimension(1:), intent(INOUT) :: header
    character (LEN=*),                 intent(IN)    :: card
    integer(kind=I4B),                 intent(IN)    :: idx
    logical(LGT),                      intent(IN)    :: overwrite
    logical(LGT),     optional,        intent(IN)    :: long_strn

    integer(i4b) :: i2, j_low, j_hi, iw, hdtype, k, step0, step
    integer(i4b) :: lencard, nlh, lenkv
    character(len=80) :: tmpline
    character (LEN=10)  :: pad10 = '    '
    logical(LGT)        :: do_long_strn
    logical(LGT), parameter :: WAC = .false. ! wrap-around long comments
    !=====================================================================
    tmpline=''
    do_long_strn = .false.
    if (present(long_strn)) do_long_strn = long_strn

    nlh = size(header)
    lencard = len_trim(card) ! total length of input card "kwd = value / comment"
    if (do_long_strn) then
       lenkv = index(card,"' / ",back=.true.) ! length of "kwd = 'value'"
       step0 = 75
    else
       lenkv = lencard
       step0 = 80
    endif

    step = step0
    j_low = 1  ! character along the card
    j_hi  = step
    k = 1
    iw = idx ! line down the header


    do ! deal with long cards
       if (j_low > lencard) exit

       if (iw > nlh) then
          print*,'WARNING: Header is too short, card not written'
          print*,trim(card)
          return
       endif

       ! leave some room for new line, unless one wants to overwrite
       if (.not.overwrite .or. k>1) then
          do i2 = nlh, iw+1, -1
             header(i2) = header(i2-1)
          enddo
       endif

       ! write line
       if (do_long_strn) then
          if (k == 1) then
             header(iw) = trim(card(j_low:j_hi))
          else
             if (j_low <= lenkv) then
                ! (long) keyword substring
                if (j_hi >= lenkv .and. WAC) then
                   header(iw) = "CONTINUE  '"//&
                        & trim(card(j_low:lenkv-1))//&
                        & "&'"//trim(card(lenkv+1:j_hi))
                else
                   header(iw) = "CONTINUE  '"//trim(card(j_low:j_hi))
                endif
             else
                ! (long) comment substring
                if (WAC) header(iw) = "CONTINUE  '&' /"//trim(card(j_low:j_hi))
             endif
          endif
          if (j_hi < lenkv) header(iw) = trim(header(iw))//"&'"
          if (j_hi < lenkv) then
             step = step0 - 10
          else
             step = step0 - 15
          endif
       else
          hdtype = 0
          statusfits = 0
          if (k == 1) then
             CALL ftgthd(card(j_low:j_hi), tmpline, hdtype, statusfits)
             header(iw) = tmpline
          else
             if (WAC) then ! wrap around long comment
                CALL ftgthd(pad10//card(j_low:j_hi), tmpline, hdtype, statusfits)
                header(iw) = tmpline
             endif
          endif
       endif

       ! select next characters from card to put in header
       j_low = j_hi + 1
       j_hi  = min(j_low + step, lencard)

       ! next header line
       k = k + 1
       iw = iw + 1

    enddo

  end subroutine insert_line_in_header
  !=====================================================================
  subroutine merge_headers( header1, header2)
    IMPLICIT none
    character (LEN=80), dimension(1:), intent(IN)    :: header1
    character (LEN=80), dimension(1:), intent(INOUT) :: header2
    integer(kind=I4B) :: iw1, iw2, s1, s2, ss

    s2 = size(header2)
    iw2 = s2
    do while(header2(iw2) == '' .and. iw2 > 1)
       iw2 = iw2 - 1
    enddo
    iw2 = iw2 + 1

    s1 = size(header1)
    iw1 = s1
    do while(header1(iw1) == '' .and. iw1 > 1)
       iw1 = iw1 - 1
    enddo
    iw1 = iw1 + 1

    ss = MIN(iw1, s2-iw2+1)

    if (ss < iw1) then
       print*,' Second header in merge_headers is not long enough'
       print*,' Should be ',iw1+iw2-2,' instead of ',s2
       print*,' Output will be truncated'
    endif

    header2(iw2:iw2+ss-1) = header1(1:ss)

    RETURN
  END subroutine merge_headers



  !===================================================================
  subroutine d_get_card(header, kwd, value, comment, count) ! double precision
    implicit none
    character (LEN=80), dimension(:), intent(IN)  :: header
    character (LEN=*),                intent(IN)  :: kwd
    real(kind=DP),                    intent(OUT) :: value
    character (LEN=*),                intent(OUT), OPTIONAL :: comment
    integer(kind=I4B),                intent(OUT), OPTIONAL :: count
    integer :: i

    count_in = 0
    value = 0.0d0
    nlh = size(header)
    do i=1, nlh ! scan header for keyword
       card = header(i)
       call ftcmps(kwd,card(1:8), casesen, match, exact)
       if (match) then
          !                        extract value as a string
          call ftpsvc(card, stval, stcom, statusfits)
          call ftdtyp(stval, dtype, statusfits) ! find its type
          if (dtype == 'F' .or. dtype == 'I') then    ! if right type
             read(stval,*) value     ! convert string to numerical value
             count_in = 1
             if (present(comment)) comment = stcom
             if (present(count))   count = count_in
             return
          else
             print*,'Uncompatible type for keyword: '//card(1:30)
             print*,'expected DOUBLE (F), found: '//dtype
             call fatal_error
          endif
       endif
    enddo
    if (verbose) print*,' >>>>> keyword '//trim(kwd)//' not found <<<<< '
    if (present(comment)) comment = ' '
    if (present(count))   count = count_in
    return
  end subroutine d_get_card

  !===================================================================
  subroutine f_get_card(header, kwd, value, comment, count) ! single precision
    implicit none
    character (LEN=80), dimension(:), intent(IN)  :: header
    character (LEN=*),                intent(IN)  :: kwd
    REAL(kind=SP),                    intent(OUT) :: value
    character (LEN=*),         intent(OUT), OPTIONAL :: comment
    integer(kind=I4B),         intent(OUT), OPTIONAL :: count
    integer :: i

    count_in = 0
    value = 0.
    nlh = size(header)
    do i=1, nlh ! scan header for keyword
       card = header(i)
       call ftcmps(kwd,card(1:8), casesen, match, exact)
       if (match) then
          !                        extract value as a string
          call ftpsvc(card, stval, stcom, statusfits)
          call ftdtyp(stval, dtype, statusfits) ! find its type
          if (dtype == 'F' .or. dtype == 'I') then    ! if right type
             read(stval,*) value     ! convert string to numerical value
             if (present(comment)) comment = stcom
             count_in = 1
             if (present(count))   count = count_in
             return
          else
             print*,'Uncompatible type for keyword: '//card(1:30)
             print*,'expected REAL (F), found: '//dtype
             call fatal_error
          endif
       endif
    enddo
    if (verbose) print*,' >>>>> keyword '//trim(kwd)//' not found <<<<< '
    if (present(comment)) comment = ' '
    if (present(count))   count = count_in
    return
  end subroutine f_get_card


  !===================================================================
  function trim_quotes(card) result(c2)
    ! remove first and last quotes
    implicit none
    character(len=*), intent(IN) :: card
    character(len=FILENAMELEN)   :: c2
    !
    character(len=FILENAMELEN)   :: c1
    integer :: ifor, ibac

    c1 = trim(adjustl(card))
    ifor = index(c1,"'")
    ibac = index(c1,"'",back=.true.)
    if (ifor >= 1                        ) c1(ifor:ifor) = " "
    if (ibac <= len(c1) .and. ibac > ifor) c1(ibac:ibac) = " "
    c2 = trim(adjustl(c1))

    return
  end function trim_quotes
  !===================================================================
  function trim_taa(card) result(c2)
    ! remove trailing Ampersand (&)
    implicit none
    character(len=*), intent(IN) :: card
    character(len=FILENAMELEN)   :: c2
    !
    character(len=FILENAMELEN)   :: c1
    integer :: ibac

    c1 = trim(adjustl(card))
    ibac = index(c1,"&",back=.true.)
    if (ibac <= len(c1)) c1(ibac:ibac) = " "
    c2 = trim(adjustl(c1))

    return
  end function trim_taa

  !===================================================================
  subroutine a_get_card(header, kwd, value, comment, count) ! ascii string
    implicit none
    character (LEN=80), dimension(:), intent(IN)  :: header
    character (LEN=*),                intent(IN) :: kwd
    character (LEN=*),                intent(OUT) :: value
    character (LEN=*),         intent(OUT), OPTIONAL :: comment
    integer(kind=I4B),         intent(OUT), OPTIONAL :: count

    character(len=80) :: stval1, stcom1, card1
    integer :: ifor, ibac, i, kaa, kq1, kq2, ip1
    logical(LGT) :: long_strn

    count_in = 0
    value = ' '
    stcom = ''
    stval1= ''
    stcom1= ''
    nlh = size(header)
    statusfits = 0
    long_strn = .false.

    ! test for possible long string keywords
    do i=1, nlh 
       card = header(i)
       call ftcmps('CONTINUE',card(1:8), casesen, long_strn, exact)
       if (long_strn) exit
    enddo


    do i=1, nlh ! scan header for keyword
       card = header(i)
       call ftcmps(kwd,card(1:8), casesen, match, exact)
       if (match) then
          if (long_strn) then ! test if current keyword is long string
             long_strn = long_strn .and. (index(card,"&'",back=.true.) > 0)
             call ftpsvc(card, stval1, stcom1, statusfits)
             value = trim(trim_taa(trim_quotes(stval1)))
             stcom = trim(stcom1)
          endif
          if (long_strn) then
             ip1 = i
             do 
                ip1 = ip1+1
                card1 = header(ip1)
                if (card1(1:10) == 'CONTINUE  ') then
                   card1(1:10)   = 'JUNK    = ' ! trick to get ftpsvc to work
                   call ftpsvc(card1, stval1, stcom1, statusfits)
                   value = trim(value)//trim(trim_taa(trim_quotes(stval1)))
                   stcom = trim(stcom)//trim(stcom1)
                else
                   exit
                endif
             enddo
          else
             !                        extract value as a string
             call ftpsvc(card, stval1, stcom1, statusfits)
             value = trim(trim_quotes(stval1))
             stcom = trim(stcom1)
          endif
          count_in = 1
          if (present(comment)) comment = stcom
          if (present(count))   count = count_in
          return
       endif
    enddo
    if (verbose) print*,' >>>>> keyword '//trim(kwd)//' not found <<<<< '
    if (present(comment)) comment = ' '
    if (present(count))   count = count_in
    return
  end subroutine a_get_card

  !===================================================================
  subroutine i_get_card(header, kwd, value, comment, count) ! integer (4Bytes)
    implicit none
    character (LEN=80), dimension(:), intent(IN)  :: header
    character (LEN=*),                intent(IN) :: kwd
    integer(kind=I4B),                          intent(OUT) :: value
    character (LEN=*),         intent(OUT), OPTIONAL :: comment
    integer(kind=I4B),         intent(OUT), OPTIONAL :: count
    integer(kind=I4B) :: i

    count_in = 0
    value = 0
    nlh = size(header)
    do i=1, nlh ! scan header for keyword
       card = header(i)
       call ftcmps(kwd,card(1:8), casesen, match, exact)
       if (match) then
          !                        extract value as a string
          call ftpsvc(card, stval, stcom, statusfits)
          call ftdtyp(stval, dtype, statusfits) ! find its type
          if (dtype == 'I') then    ! if right type
             read(stval,*) value     ! convert string to numerical value
             if (present(comment)) comment = stcom
             count_in = 1
             if (present(count))   count = count_in
             return
          else
             print*,'Uncompatible type for keyword: '//card(1:30)
             print*,'expected integer (I), found: '//dtype
             call fatal_error
          endif
       endif
    enddo
    if (verbose) print*,' >>>>> keyword '//trim(kwd)//' not found <<<<< '
    if (present(comment)) comment = ' '
    if (present(count))   count = count_in
    return
  end subroutine i_get_card

  !===================================================================
#ifndef NO64BITS
  ! compile 8 byte integer flavor only if possible
  subroutine ii_get_card(header, kwd, value, comment, count) ! integer (8Bytes)
    implicit none
    character (LEN=80), dimension(:), intent(IN)  :: header
    character (LEN=*),                intent(IN) :: kwd
    integer(kind=I8B),                          intent(OUT) :: value
    character (LEN=*),         intent(OUT), OPTIONAL :: comment
    integer(kind=I4B),         intent(OUT), OPTIONAL :: count
    integer(kind=I4B) :: i

    count_in = 0
    value = 0
    nlh = size(header)
    do i=1, nlh ! scan header for keyword
       card = header(i)
       call ftcmps(kwd,card(1:8), casesen, match, exact)
       if (match) then
          !                        extract value as a string
          call ftpsvc(card, stval, stcom, statusfits)
          call ftdtyp(stval, dtype, statusfits) ! find its type
          if (dtype == 'I') then    ! if right type
             read(stval,*) value     ! convert string to numerical value
             if (present(comment)) comment = stcom
             count_in = 1
             if (present(count))   count = count_in
             return
          else
             print*,'Uncompatible type for keyword: '//card(1:30)
             print*,'expected integer (I), found: '//dtype
             call fatal_error
          endif
       endif
    enddo
    if (verbose) print*,' >>>>> keyword '//trim(kwd)//' not found <<<<< '
    if (present(comment)) comment = ' '
    if (present(count))   count = count_in
    return
  end subroutine ii_get_card
#endif
  !===================================================================
  subroutine l_get_card(header, kwd, value, comment, count) ! logical
    implicit none
    character (LEN=80), dimension(:), intent(IN)  :: header
    character (LEN=*),                intent(IN)  :: kwd
    logical(kind=LGT),                intent(OUT) :: value
    character (LEN=*),         intent(OUT), OPTIONAL :: comment
    integer(kind=I4B),         intent(OUT), OPTIONAL :: count
    integer(kind=I4B) :: i

    count_in = 0
    value = .false.
    nlh = size(header)
    do i=1, nlh ! scan header for keyword
       card = header(i)
       call ftcmps(kwd,card(1:8), casesen, match, exact)
       if (match) then
          !                        extract value as a string
          call ftpsvc(card, stval, stcom, statusfits)
          call ftdtyp(stval, dtype, statusfits) ! find its type
          if (dtype == 'L') then    ! if right type
             read(stval,*) value     ! convert string to numerical value
             if (present(comment)) comment = stcom
             count_in = 1
             if (present(count))   count = count_in
             return
          else
             print*,'Uncompatible type for keyword: '//card(1:30)
             print*,'expected logical (L), found: '//dtype
             call fatal_error
          endif
       endif
    enddo
    if (verbose) print*,' >>>>> keyword '//trim(kwd)//' not found <<<<< '
    if (present(comment)) comment = ' '
    if (present(count))   count = count_in
    return
  end subroutine l_get_card

  !===================================================================
  subroutine del_cardn(header, kwds) ! remove cards
    !===================================================================
    ! remove for header the card corresponding to keyword kwds
    ! kwds can be a vector
    !===================================================================
    implicit none
    character(len=80), dimension(1:), intent(INOUT) :: header
    character(len=*),  dimension(1:), intent(IN)    :: kwds
    integer(kind=I4B) :: i
    character(len=20) :: kwd
    !===================================================================
    nlh = size(kwds)
    do i = 1, nlh
       kwd = adjustl(kwds(i))
       if (trim(kwd) /= "") then
          kwd = "- "//kwd
          call write_hl(header, kwd)
       endif
    enddo
  end subroutine del_cardn
  !===================================================================
  subroutine del_card1(header, kwds) ! remove cards
    !===================================================================
    ! remove for header the card corresponding to keyword kwds
    ! kwds can be a vector
    !===================================================================
    implicit none
    character(len=80), dimension(1:), intent(INOUT) :: header
    character(len=*),                 intent(IN)    :: kwds
    character(len=20) :: kwd
    !===================================================================
       kwd = adjustl(kwds)
       if (trim(kwd) /= "") then
          kwd = "- "//kwd
          call write_hl(header, kwd)
       endif
  end subroutine del_card1


  !===================================================================
  subroutine write_minimal_header(header &
       &, dtype & ! 'alm', 'cl', 'cutmap', 'map'
       &, append &
       &, nside &            ! required for 'map', 'cutmap'
       &, order, ordering &  ! required for 'map', 'cutmap'
       &, coordsys &         ! required for 'map', 'cutmap'
       &, creator &
       &, version &
       &, randseed &
       &, beam_leg &
       &, fwhm_degree &
       &, units &
       &, nlmax &
       &, polar &
       &, nmmax &            ! relevant for 'alm'
       &, bcross &           ! relevant for 'cl'
       &, deriv &            ! relevant for 'map'
       &, asym_cl &          ! relevant for 'cl'
)
  !===================================================================
    
    implicit none
    character(len=*), intent(inout), dimension(:) :: header
    character(len=*), intent(in)           :: dtype
    integer(i4b),     intent(in), optional :: order, randseed, deriv
!    integer(i4b), intent(in), optional :: extension
    integer(i4b),     intent(in), optional :: nside, nlmax, nmmax
    logical(LGT),     intent(in), optional :: append, polar, bcross, asym_cl
    character(len=*), intent(in), optional :: creator, version, coordsys, ordering
    character(len=*), intent(in), optional :: beam_leg
!     character(len=*), intent(in), dimension(1:), optional :: units
    character(len=*), intent(in), optional :: units
    real(dp),         intent(in), optional :: fwhm_degree

    integer(i4b) :: ext, it, iq, iu, my_deriv
    character(len=80) :: my_coordsys
    character(len=80), dimension(1:20) :: my_units
    character(len=8) :: my_ordering
    character(len=FILENAMELEN) :: udtype, my_beam_leg
    logical(LGT) :: do_full, do_alm, do_cl, do_cut, do_partial, &
         & do_polar, do_bcross, do_reset, do_asym, do_map

    character(len=*), parameter :: code = 'write_minimal_header'
    !--------------------------------------------------------------------------------
    ! find out type of data to be contained in FITS file
    udtype = strupcase(dtype)
    do_alm = .false. ; do_cl  = .false. ; do_full = .false. 
    do_cut = .false. ; do_partial = .false.
    select case (trim(udtype))
    case ('ALM')
       do_alm = .true.
    case ('CL')
       do_cl = .true.
    case ('MAP')
       do_full = .true.
    case ('CUTMAP')
       do_cut = .true.
    ! case ('PARTIAL')
       do_partial = .true.
    case default
       print*,'Invalid choice of datatype:'//trim(dtype)
       print*,'Should be one of: ALM, alm, CL, cl, CUTMAP, cutmap, MAP, map'
       call fatal_error(code)
    end select
    do_map = (do_cut .or. do_full .or. do_partial)

    !
    if (present(order) .and. present(ordering)) then
       print*,'Choose either ORDER (1/2) or ORDERING (RING/NESTED)'
       call fatal_error(code)
    endif

    ! check for required keywords for map
    if (do_map) then
       if (.not. present(nside)) then
          print*,'NSIDE required for maps in '//code
          call fatal_error(code)
       endif
       if (.not.present(order) .and. .not.present(ordering)) then
          print*,'ORDERING required for maps in '//code
          call fatal_error(code)
       endif
    endif

       
    ! set default values
    do_reset = .true.
    if (present(append)) do_reset = .not. append
    my_coordsys = 'unknown'
    if (present(coordsys)) then
       if (trim(coordsys) /= '') my_coordsys = coordsys
    endif
    my_units(:) = 'unknown'
    if (present(order)) then
       my_ordering = "RING"
       if (order /=1 .and. order /=2) then
          print*,'ORDER should be 1 or 2, it is: ',order
          call fatal_error(code)
       endif
       if (order == 2)  my_ordering = "NESTED"
    endif
    if (present(ordering)) then
       my_ordering = strupcase(ordering)
       if (trim(my_ordering) /= 'RING' .and. trim(my_ordering) /= 'NESTED') then
          print*,'ORDERING should be RING or NESTED, it is: ',trim(ordering)
          call fatal_error(code)
       endif
    endif
    do_bcross = .false.
    if (present(bcross)) then
       do_bcross = bcross
    endif
    do_asym = .false.
    if (present(asym_cl)) then
       do_asym = asym_cl
    endif
    my_deriv = 0
    if (present(deriv)) my_deriv = deriv
    if (my_deriv <0 .or. my_deriv >2) then
       print*,'DERIV out of range [0,2]: ',my_deriv
       call fatal_error(code)
    endif
    my_beam_leg = ''
    if (present(beam_leg)) my_beam_leg = beam_leg

    ! reset header if requested
    if (do_reset) header(:) = " "

    do_polar = .false.
    if (present(polar)) do_polar = polar
    
    ext = 0
!     if (present(extension)) ext = extension

    if (do_alm) then
       !--------------------------------------------------------------------------------------------
       !                            minimal header for ALM
       !--------------------------------------------------------------------------------------------
       ! missing : CUT-SKY, NITERALM, LREGRESS
!        select case (ext)
!        case (0)
!           call add_card(header,"EXTNAME","'a_lms (TEMPERATURE)'")
!        case (1)
!           call add_card(header,"EXTNAME","'a_lms (GRAD / ELECTRIC component)'")
!        case (2)
!           call add_card(header,"EXTNAME","'a_lms (CURL / MAGNETIC component)'")
!        end select
       call add_card(header,"COMMENT","-----------------------------------------------")
       call add_card(header,"COMMENT","     Map Analysis Specific Keywords      ")
       call add_card(header,"COMMENT","-----------------------------------------------")
       call add_card(header,"EXTNAME","'a_lms coefficients'") ! place holder
       call add_card(header)
       if (present(creator)) call add_card(header,"CREATOR",creator, "Software creating the FITS file")
       if (present(version)) call add_card(header,"VERSION",version, "Version of the simulation software")
       call add_card(header)
       if (present(nlmax)) call add_card(header,"MAX-LPOL",nlmax  ,"Maximum Legendre order L")
       if (present(nmmax)) call add_card(header,"MAX-MPOL",nmmax  ,"Maximum Legendre degree M")
       if (present(randseed))  call add_card(header,"RANDSEED",randseed,"Random generator seed")
       if (trim(my_beam_leg)/='') then
          call add_card(header,"BEAM_LEG",trim(my_beam_leg),"File containing Legendre transform of symmetric beam")
       else if (present(fwhm_degree)) then
          call add_card(header,"FWHM"    ,fwhm_degree   ," [deg] FWHM of gaussian symmetric beam")
       endif
       call add_card(header)

       
       if (present(units)) then
!           iu = min(ext+1, size(units))
!           my_units(1) = "'" // trim(units(iu)) // "'"
          my_units(1) = "'" // trim(units) // "'"
       endif
       call add_card(header,"COMMENT"," The real and imaginary part of the a_lm with m>=0")
       call add_card(header)
       call add_card(header,"TTYPE1", "INDEX"," i = l^2 + l + m + 1")
       call add_card(header,"TUNIT1", "   "," index")
       call add_card(header)
       call add_card(header,"TTYPE2", "REAL"," REAL a_lm")
       call add_card(header,"TUNIT2", my_units(1)," alm units")
       call add_card(header)
       !
       call add_card(header,"TTYPE3", "IMAG"," IMAGINARY a_lm")
       call add_card(header,"TUNIT3", my_units(1)," alm units")
       call add_card(header)

       if (present(nside)) call add_card(header,"NSIDE"   ,nside,   "Resolution parameter for HEALPIX")
       call add_card(header,"COORDSYS",trim(my_coordsys),"Input map coordinate system")
       call add_card(header,"COMMENT","G = Galactic, E = ecliptic, C = celestial = equatorial")
       call add_card(header,"POLAR",do_polar," Polarisation included (True/False)")
    endif

    if (do_cl) then
       !--------------------------------------------------------------------------------------------
       !                            minimal header for CL
       !--------------------------------------------------------------------------------------------
       ! missing: CUT-SKY, NITERALM, LREGRESS, NMAPS_IN
       call add_card(header,"COMMENT","-----------------------------------------------")
       call add_card(header,"COMMENT","     Map Analysis Specific Keywords      ")
       call add_card(header,"COMMENT","-----------------------------------------------")
       call add_card(header,"EXTNAME","'POWER SPECTRUM'")
       call add_card(header)
       if (present(creator)) call add_card(header,"CREATOR",creator, "Software creating the FITS file")
       if (present(version)) call add_card(header,"VERSION",version, "Version of the simulation software")
       call add_card(header)
       if (present(nlmax)) call add_card(header,"MAX-LPOL",nlmax  ,"Maximum Legendre order L")
       if (present(nmmax)) call add_card(header,"MAX-MPOL",nmmax  ,"Maximum Legendre degree M")
       if (present(randseed))  call add_card(header,"RANDSEED",randseed,"Random generator seed")
       call add_card(header)
       call add_card(header,"POLAR",do_polar," Polarisation included (True/False)")
       call add_card(header,"BCROSS",do_bcross," Magnetic cross terms included (True/False)")
       call add_card(header,"ASYMCL",do_asym," Asymmetric pol cross terms (XY vs YX) included")
       call add_card(header)
       if (present(nside)) call add_card(header,"NSIDE",nside," Resolution Parameter of Input Map")
       call add_card(header,"COORDSYS",my_coordsys," Input Map Coordinate system")
       call add_card(header,"COMMENT","G = Galactic, E = ecliptic, C = celestial = equatorial")
       call add_card(header)
       call add_card(header,"TTYPE1", "TEMPERATURE","Temperature C(l)")
!        if (present(units)) my_units(1) = "'"//trim(units(1))//"'" ! add quotes to keep words together
       if (present(units)) my_units(1) = "'"//trim(units)//"'" ! add quotes to keep words together
       call add_card(header,"TUNIT1", my_units(1),"square of map units")
       call add_card(header)
       !
       if (do_polar) then
          call add_card(header,"TTYPE2", "GRADIENT","Gradient (=ELECTRIC) polarisation C(l)")
          call add_card(header,"TUNIT2", my_units(1),"power spectrum units")
          call add_card(header)
          !
          call add_card(header,"TTYPE3", "CURL","Curl (=MAGNETIC) polarisation C(l)")
          call add_card(header,"TUNIT3", my_units(1),"power spectrum units")
          call add_card(header)
          !
          !call add_card(header,"TTYPE4", "G-T","Gradient-Temperature (=CROSS) terms")
          call add_card(header,"TTYPE4", "TG","Temperature X Gradient cross terms")
          call add_card(header,"TUNIT4", my_units(1),"power spectrum units")
          call add_card(header)
          !
          if (do_bcross) then
             !call add_card(header,"TTYPE5", "C-T","Curl-Temperature terms")
             call add_card(header,"TTYPE5", "TC","Temperature X Curl cross terms")
             call add_card(header,"TUNIT5", my_units(1),"power spectrum units")
             call add_card(header)
           !
             !call add_card(header,"TTYPE6", "C-G","Curl-Gradient terms")
             call add_card(header,"TTYPE6", "GC","Gradient x Curl cross terms")
             call add_card(header,"TUNIT6", my_units(1),"power spectrum units")
             call add_card(header)
           !
             if (do_asym) then
                call add_card(header,"TTYPE7", "GT","Gradient1 x Temperature2 cross terms")
                call add_card(header,"TUNIT7", my_units(1),"power spectrum units")
                call add_card(header)

                call add_card(header,"TTYPE8", "CT","Curl1 x Temperature2 cross terms")
                call add_card(header,"TUNIT8", my_units(1),"power spectrum units")
                call add_card(header)

                call add_card(header,"TTYPE9", "CG","Curl1 x Gradient2 cross terms")
                call add_card(header,"TUNIT9", my_units(1),"power spectrum units")
                call add_card(header)
             endif
          else
             if (do_asym) then
                call add_card(header,"TTYPE5", "GT","Gradient1 x Temperature2 cross terms")
                call add_card(header,"TUNIT5", my_units(1),"power spectrum units")
                call add_card(header)
             endif
          endif
          !
          call add_card(header,"COMMENT","The polarisation power spectra have the same definition as in CMBFAST")
          call add_card(header)
       endif

    endif ! cl

    if (do_map) then
       !--------------------------------------------------------------------------------------------
       !                            minimal header for FULL SKY MAP and CUT SKY MAP
       !--------------------------------------------------------------------------------------------
       call add_card(header)
       call add_card(header,"COMMENT","-----------------------------------------------")
       call add_card(header,"COMMENT","     Sky Map Pixelisation Specific Keywords    ")
       call add_card(header,"COMMENT","-----------------------------------------------")
       call add_card(header,"PIXTYPE","HEALPIX","HEALPIX Pixelisation")
       call add_card(header,"ORDERING",my_ordering,  "Pixel ordering scheme, either RING or NESTED")
       
       call add_card(header,"NSIDE"   ,nside,   "Resolution parameter for HEALPIX")
       if (do_full) then
          call add_card(header,"FIRSTPIX",0,"First pixel # (0 based)")
          call add_card(header,"LASTPIX",(12_i8b*nside)*nside-1,"Last pixel # (0 based)")
       endif
       call add_card(header,"COORDSYS",trim(my_coordsys),"Pixelisation coordinate system")
       call add_card(header,"COMMENT","G = Galactic, E = ecliptic, C = celestial = equatorial")
       call add_card(header,"BAD_DATA",  HPX_DBADVAL ,"Sentinel value given to bad pixels")
       call add_card(header) ! blank line


       call add_card(header,"COMMENT","-----------------------------------------------")
       call add_card(header,"COMMENT","     Planck Simulation Specific Keywords      ")
       call add_card(header,"COMMENT","-----------------------------------------------")
       call add_card(header,"EXTNAME","MAP") ! place holder, to be updated
       if (present(creator)) call add_card(header,"CREATOR",creator, "Software creating the FITS file")
       if (present(version)) call add_card(header,"VERSION",version, "Version of the simulation software")
       call add_card(header)
       if (present(nlmax)) call add_card(header,"MAX-LPOL",nlmax  ,"Maximum Legendre order L")
       if (present(nmmax)) call add_card(header,"MAX-MPOL",nmmax  ,"Maximum Legendre degree M")
       if (present(randseed))  call add_card(header,"RANDSEED",randseed,"Random generator seed")
       call add_card(header,"POLCCONV","COSMO","Coord. convention for polarisation (COSMO/IAU)")
       if (trim(my_beam_leg)/='') then
          call add_card(header,"BEAM_LEG",trim(my_beam_leg),"File containing Legendre transform of symmetric beam")
       else if (present(fwhm_degree)) then
          call add_card(header,"FWHM"    ,fwhm_degree   ," [deg] FWHM of gaussian symmetric beam")
       endif

       call add_card(header,"COMMENT","-----------------------------------------------")
       call add_card(header,"COMMENT","     Data Description Specific Keywords       ")
       call add_card(header,"COMMENT","-----------------------------------------------")
       call add_card(header)
       if (do_full) then
          call add_card(header,"EXTNAME","'FULL SKY MAP'",update=.true.)
          call add_card(header,"COMMENT","Full sky data")
          call add_card(header,"OBJECT","FULLSKY")
          call add_card(header,"INDXSCHM","IMPLICIT"," Indexing : IMPLICIT or EXPLICIT")
!          call add_card(header,"GRAIN", 0, " Grain of pixel indexing")
       endif
       if (do_cut .or. do_partial)  then
          call add_card(header,"EXTNAME","'CUT SKY MAP'",update=.true.)
          call add_card(header,"COMMENT","Cut sky data")
          call add_card(header,"OBJECT","PARTIAL")
          call add_card(header,"INDXSCHM","EXPLICIT"," Indexing : IMPLICIT or EXPLICIT")
!          call add_card(header,"GRAIN", 1, " Grain of pixel indexing")
       endif
!        call add_card(header,"COMMENT","GRAIN=0 : no indexing of pixel data                         (IMPLICIT)")
!        call add_card(header,"COMMENT","GRAIN=1 : 1 pixel index -> 1 pixel data                     (EXPLICIT)")
!        call add_card(header,"COMMENT","GRAIN>1 : 1 pixel index -> data of GRAIN consecutive pixels (EXPLICIT)")
       call add_card(header) ! blank line
       call add_card(header,"POLAR",do_polar," Polarisation included (True/False)")
       if (do_cut) then
          call add_card(header) ! blank line
          call add_card(header,"TTYPE1", "PIXEL","pixel index")
          call add_card(header,"TUNIT1", "",     "")
          call add_card(header)
          if (present(units)) then
             my_units = units
             call add_card(header,"TUNIT2", my_units(1),"physical units of map")
             call add_card(header,"TUNIT4", my_units(1),"physical units of map")
          endif
       endif
!        if (do_partial) then
!           if (present(units)) my_units = units
!           it = 1 ; iq = 1 ; iu = 1
!           call add_card(header) ! blank line
!           call add_card(header,"TTYPE1", "PIXEL","pixel index")
!           call add_card(header,"TUNIT1", "",     "")
!           call add_card(header)

!           call add_card(header) ! blank line
!           call add_card(header,"TTYPE2", "TEMPERATURE","Temperature map")
!           call add_card(header,"TUNIT2", my_units(it), "physical units of map")
!           call add_card(header)

!           if (do_polar) then
!              call add_card(header,"TTYPE3", "Q_POLARISATION","Q Polarisation map")
!              call add_card(header,"TUNIT3", my_units(iq),    "physical units of map")
!              call add_card(header)

!              call add_card(header,"TTYPE4", "U_POLARISATION","U Polarisation map")
!              call add_card(header,"TUNIT4", my_units(iu),    "physical units unit")
!              call add_card(header)
!           endif
!        endif
       if (do_full) then
          call add_card(header,"DERIV",my_deriv," Derivative included (0, 1 or 2)")

          if (present(units)) then
!              do it=1, size(units)
!                 my_units(it) = "'" // trim(units(it)) // "'"
!              enddo
!              it = min(1, size(units))
!              iq = min(2, size(units))
!              iu = min(3, size(units))
             my_units = units
             it = 1 ; iq = 1 ; iu = 1
          else
             it = 1 ; iq = 1 ; iu = 1
          endif
          
          call add_card(header) ! blank line
          call add_card(header,"TTYPE1", "TEMPERATURE","Temperature map")
          call add_card(header,"TUNIT1", my_units(it),"map unit")
          call add_card(header)

          if (do_polar) then
             call add_card(header,"TTYPE2", "Q_POLARISATION","Q Polarisation map")
             call add_card(header,"TUNIT2", my_units(iq),"map unit")
             call add_card(header)

             call add_card(header,"TTYPE3", "U_POLARISATION","U Polarisation map")
             call add_card(header,"TUNIT3", my_units(iu),"map unit")
             call add_card(header)

             if (my_deriv>=1) then
                call add_card(header,"TTYPE4", "DX_T","1st derivative of T wrt theta")
                call add_card(header,"TUNIT4", my_units(it),"map unit")
                call add_card(header)
                call add_card(header,"TTYPE5", "DX_Q","1st derivative of Q wrt theta")
                call add_card(header,"TUNIT5", my_units(iq),"map unit")
                call add_card(header)
                call add_card(header,"TTYPE6", "DX_U","1st derivative of U wrt theta")
                call add_card(header,"TUNIT6", my_units(iu),"map unit")
                call add_card(header)

                call add_card(header,"TTYPE7", "DY_T","(1st derivative of T wrt phi)/sin(theta)")
                call add_card(header,"TUNIT7", my_units(it),"map unit")
                call add_card(header)
                call add_card(header,"TTYPE8", "DY_Q","(1st derivative of Q wrt phi)/sin(theta)")
                call add_card(header,"TUNIT8", my_units(iq),"map unit")
                call add_card(header)
                call add_card(header,"TTYPE9", "DY_U","(1st derivative of U wrt phi)/sin(theta)")
                call add_card(header,"TUNIT9", my_units(iu),"map unit")
                call add_card(header)
             endif

             if (my_deriv>=2) then
                call add_card(header,"TTYPE10", "DXX_T","2nd derivative of T wrt theta, theta")
                call add_card(header,"TUNIT10", my_units(it),"map unit")
                call add_card(header)
                call add_card(header,"TTYPE11", "DXX_Q","2nd derivative of Q wrt theta, theta")
                call add_card(header,"TUNIT11", my_units(iq),"map unit")
                call add_card(header)
                call add_card(header,"TTYPE12", "DXX_U","2nd derivative of U wrt theta, theta")
                call add_card(header,"TUNIT12", my_units(iu),"map unit")
                call add_card(header)

                call add_card(header,"TTYPE13", "DXY_T","(2nd deriv of T wrt theta, phi)/sin(theta)")
                call add_card(header,"TUNIT13", my_units(it),"map unit")
                call add_card(header)
                call add_card(header,"TTYPE14", "DXY_Q","(2nd deriv of Q wrt theta, phi)/sin(theta)")
                call add_card(header,"TUNIT14", my_units(iq),"map unit")
                call add_card(header)
                call add_card(header,"TTYPE15", "DXY_U","(2nd deriv of U wrt theta, phi)/sin(theta)")
                call add_card(header,"TUNIT15", my_units(iu),"map unit")
                call add_card(header)

                call add_card(header,"TTYPE16", "DYY_T","(2nd deriv of T wrt phi, phi)/sin(theta)^2")
                call add_card(header,"TUNIT16", my_units(it),"map unit")
                call add_card(header)
                call add_card(header,"TTYPE17", "DYY_Q","(2nd deriv of Q wrt phi, phi)/sin(theta)^2")
                call add_card(header,"TUNIT17", my_units(iq),"map unit")
                call add_card(header)
                call add_card(header,"TTYPE18", "DYY_U","(2nd deriv of U wrt phi, phi)/sin(theta)^2")
                call add_card(header,"TUNIT18", my_units(iu),"map unit")
                call add_card(header)
             endif
          else ! not polarized
             if (my_deriv>=1) then
                call add_card(header,"TTYPE2", "DX_T","1st derivative of T wrt theta")
                call add_card(header,"TUNIT2", my_units(it),"map unit")
                call add_card(header)

                call add_card(header,"TTYPE3", "DY_T","(1st derivative of T wrt phi)/sin(theta)") ! corrected 2009-01-08
                call add_card(header,"TUNIT3", my_units(it),"map unit")
                call add_card(header)
             endif

             if (my_deriv>=2) then
                call add_card(header,"TTYPE4", "DXX_T","2nd derivative of T wrt theta, theta")
                call add_card(header,"TUNIT4", my_units(it),"map unit")
                call add_card(header)

                call add_card(header,"TTYPE5", "DXY_T","(2nd deriv wrt of T theta, phi)/sin(theta)")
                call add_card(header,"TUNIT5", my_units(it),"map unit")
                call add_card(header)

                call add_card(header,"TTYPE6", "DYY_T","(2nd deriv wrt of T phi, phi)/sin(theta)^2")
                call add_card(header,"TUNIT6", my_units(it),"map unit")
                call add_card(header)
             endif
          endif ! polarized
       endif ! full sky map
    endif ! map
    

  end subroutine write_minimal_header

end module head_fits
