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
! template for routine SP/DP overloading for module fitstools
!
!    subroutine input_map_KLOAD            (4/8)
!    subroutine read_bintab_KLOAD          (4/8)
!    subroutine read_conbintab_KLOAD           NotYet
!    subroutine write_bintab_KLOAD         (4/8)
!    subroutine write_asctab_KLOAD          NA
!    subroutine dump_alms_KLOAD                NotYet
!    subroutine write_alms_KLOAD               NotYet
!    subroutine read_alms_KLOAD                NotYet
!    subroutine read_bintod_KLOAD           (8)
!    subroutine write_bintabh_KLOAD         (8)
!    subroutine unfold_weights_KLOAD
!    subroutine unfold_weightslist_KLOAD
!    subroutine read_fits_partial_KLOAD     (4/8)
!    subroutine write_fits_partial_KLOAD     (4)
!
!
! K M A P   : map kind                 either SP or DP
!
! K L O A D : suffixe of routine name, to be replaced by either s (SP) or d (DP)
!
! edited Jan 11, 2006 to deal correctly with polarised alms in write_alms and alms2fits
! edited Apr 04, 2006 to avoid concatenation problem in TFORMs  (write_plm and write_bintabh with Ifort 9)
! edited Dec 20, 2006 to accept alm file with less elements than expected (in particular if vanishing ones are not included)
! 2007-05-10 : increased maxdim (max number of columns per extension) from 20 to 40
! 2008-08-27 : in dump_alms and write_alms and write_*tab*: 
!  do not write TTYPE# and TFORM# in excess of # of fields in the file
! 2008-10-14: corrected bug introduced in write_asctab
! 2012-02-23: correction of a possible bug with index writing in dump_alms and write_alms
! 2013-12-13: increased MAXDIM from 40 to MAXDIM_TOP
! 2018-05-22: added unfold_weights, unfold_weightslist
! 2019-10-14: can write TTYPE??? (ie with up to 3 digits) in write_* and dump_alms
!

!=======================================================================
! map_bad_pixels(map, fin, fout, nbads, verbose)
!   map: input data (2D)
!   fin: input value of 'bad' pixels
!   fout: output value of same 'bad' pixels
!   nbads: number of bad pixels found
!   verbose: OPTIONAL, if set, be verbose
!=======================================================================

subroutine map_bad_pixels_KLOAD(map, fin, fout, nbads, verbose)
  use long_intrinsic, only: long_size
  real(KMAP), dimension(0:,1:), intent(inout) :: map
  real(KMAP),                   intent(in)    :: fin, fout
  integer(I8B),                 intent(out)   :: nbads
  logical(LGT),     optional,   intent(in)    :: verbose
  integer(I8B) :: i, npix
  integer(I8B), dimension(1:100) :: imiss
  integer(I4B) :: imap, nmaps
  logical(LGT) :: be_verbose
  real(KMAP) :: threshold
  !-----------------------------------------
  
  npix  = long_size(map, 1)
  nmaps = long_size(map, 2)
  threshold = abs(fin * 1.e-5_KMAP)
  
  imiss(:) = 0
  do imap = 1, nmaps
     do i=0,npix-1
        if ( abs( map(i,imap) - fin ) < threshold ) then
           map(i,imap) = fout
           imiss(imap) = imiss(imap)+1
        endif
     enddo
  enddo
  nbads = sum(imiss)
  
  be_verbose = .false.
  if (present(verbose)) be_verbose=verbose
  if (be_verbose) then
     write(*,'(a,1pe11.4)') 'blank value : ' ,fin
     do imap = 1, nmaps
        if (imiss(imap) > 0) then
           write(*,'(i12,a,f7.3,a,1pe11.4)') &
                &           imiss(imap),' missing pixels (', &
                &           (100.0_KMAP*imiss(imap))/npix,' %),'// &
                &           ' have been set to : ',fout
        endif
     enddo
  endif
  
  return
end subroutine map_bad_pixels_KLOAD
!=======================================================================
! input_map
!     reads fits file
!     filename = fits file (input)
!     map      = data read from the file (ouput) = real*4 array of size (npixtot,nmaps)
!     npixtot  = length of the map (input)
!     nmaps     = number of maps
!     fmissval  = OPTIONAL argument (input) with the value to be given to missing
!             pixels, its default value is 0
!     header    = OPTIONAL (output) contains extension header
!     units     = OPTIONAL (output) contains the maps units
!     extno     = OPTIONAL (input)  contains the unit number to read from (0 based)
!     ignore_polcconv = OPTIONAL (input), LGT, default=.false.
!                  take into account or not the POLCCONV FITS keyword
!
!     modified Feb 03 for units argument to run with Sun compiler
!     2017-02-15: deals with POLCCONV
!=======================================================================
#ifndef NO64BITS
subroutine input_map4_KLOAD(filename, map, npixtot, nmaps, &
     & fmissval, header, units, extno, ignore_polcconv)
    !=======================================================================

    INTEGER(I4B),     INTENT(IN)                           :: npixtot
    INTEGER(I4B),     INTENT(IN)                           :: nmaps
    REAL(KMAP),       INTENT(OUT), dimension(0:,1:)        :: map
    CHARACTER(LEN=*), INTENT(IN)                           :: filename
    REAL(KMAP),       INTENT(IN),                 OPTIONAL :: fmissval
    CHARACTER(LEN=*), INTENT(OUT), dimension(1:), OPTIONAL :: header
    CHARACTER(LEN=*), INTENT(OUT), dimension(1:), OPTIONAL :: units
    INTEGER(I4B),     INTENT(IN)                , optional :: extno
    logical(LGT),     intent(IN)                , optional :: ignore_polcconv
    integer(i8b) :: npixtot8

    npixtot8 = npixtot
    call input_map8_KLOAD(filename, map, npixtot8, nmaps, &
         fmissval, header, units, extno, ignore_polcconv)

    return
  end subroutine input_map4_KLOAD
#endif

!=======================================================================
subroutine input_map8_KLOAD(filename, map, npixtot, nmaps, &
     fmissval, header, units, extno, ignore_polcconv)
    !=======================================================================
  use head_fits, only: get_card, add_card
  use pix_tools, only: nside2npix
    INTEGER(I8B),     INTENT(IN)                           :: npixtot
    INTEGER(I4B),     INTENT(IN)                           :: nmaps
    REAL(KMAP),       INTENT(OUT), dimension(0:,1:)        :: map
    CHARACTER(LEN=*), INTENT(IN)                           :: filename
    REAL(KMAP),       INTENT(IN),                 OPTIONAL :: fmissval
    CHARACTER(LEN=*), INTENT(OUT), dimension(1:), OPTIONAL :: header
    CHARACTER(LEN=*), INTENT(OUT), dimension(1:), OPTIONAL :: units
    INTEGER(I4B),     INTENT(IN)                , optional :: extno
    logical(LGT),     intent(IN)                , optional :: ignore_polcconv

    INTEGER(I8B) :: i, imissing, obs_npix, maxpix, minpix
    REAL(KMAP)   :: fmissing, fmiss_effct
    integer(I4B) :: imap, nlheader

    LOGICAL(LGT) :: anynull, do_polcconv
    ! Note : read_fits_cut4 still SP and I4B only
    integer(I4B), dimension(:),   allocatable :: pixel4
    real(SP),     dimension(:),   allocatable :: signal
    ! while read_fits_partial is SP/DP and I4B/I8B
    integer(I8B), dimension(:),   allocatable :: pixel8
    real(KMAP),   dimension(:,:), allocatable :: iqu
    integer(I4B) :: status
    integer(I4B) :: type_fits, nmaps_fits
    CHARACTER(LEN=80)  :: units1
    CHARACTER(LEN=80),dimension(1:30)  :: unitsm
    character(len=80), dimension(:), allocatable  :: hbuffer
!    CHARACTER(LEN=80),dimension(:), allocatable  :: unitsm
    integer(I4B) :: extno_i, extno_f, n_ext
    CHARACTER(LEN=80)  :: polcconv
    integer(I4B)       :: ipolcconv, polar, nside_fits
    character(len=*), parameter :: primer_url = 'http://healpix.sf.net/pdf/intro.pdf'
    character(len=*), parameter :: code = 'input_map'
    
    !-----------------------------------------------------------------------

    units1    = ' '
    unitsm(:) = ' '
    fmiss_effct = 0.
    imissing = 0
    if (PRESENT(fmissval)) fmiss_effct = fmissval
    if (PRESENT(header)) then
       nlheader = size(header)
    else
       nlheader = 36*100
    endif
    extno_i = 0
    if (present(extno)) extno_i = extno
    allocate(hbuffer(1:nlheader))
    do_polcconv = .true.
    if (present(ignore_polcconv)) then
       do_polcconv = .not. ignore_polcconv
    endif

    n_ext = getnumext_fits(filename)
    obs_npix = getsize_fits(filename, nmaps = nmaps_fits, type=type_fits, extno=extno_i, &
         polcconv=ipolcconv, polarisation=polar, nside=nside_fits)

    !       if (nmaps_fits > nmaps) then 
    !          print*,trim(filename)//' only contains ',nmaps_fits,' maps'
    !          print*,' You try to read ',nmaps
    !       endif
    if (type_fits == 0 .or. type_fits == 2) then ! full sky map (in image or binary table)
       call read_bintab(filename, map(0:,1:), &
            & npixtot, nmaps, fmissing, anynull, header=hbuffer(1:), &
            & units=unitsm(1:), extno=extno_i)
       if (present(header)) then
          do i=1,nlheader
             header(i) = hbuffer(i)
          enddo
       endif
       if (present(units)) then
          units(1:size(units)) = unitsm(1:size(units))
       endif
       call map_bad_pixels(map, fmissing, fmiss_effct, imissing, verbose=.false.)
       if (imissing > 0) write(*,'(a,1pe11.4)') 'blank value : ' ,fmissing
    else if (type_fits == 3 .and. (nmaps == 1 .or. nmaps == 3) &
         & .and. (nmaps_fits > nmaps) .and. n_ext == 1) then
       ! partial FITS file, reading 1 map (I) or 3 maps (I,Q,U) from the same extension
       obs_npix = getsize_fits(filename)
       allocate(pixel8(0:obs_npix-1),          stat = status)
       call assert_alloc(status, code, "pixel8")
       allocate(iqu   (0:obs_npix-1, 1:nmaps), stat = status)
       call assert_alloc(status, code, "iqu")
       call read_fits_partial(filename, pixel8, iqu, header=hbuffer)
       if (present(header)) then
          do i=1,nlheader
             header(i) = hbuffer(i)
          enddo
       endif
       if (present(units)) then
!          do imap=1,min(size(units),nmaps-1) ! 2022-06-10
          do imap=1,min(size(units),nmaps)
             call get_card(hbuffer,'TUNIT'//string(imap+1,format='(i0.0)'), units(imap))
          enddo
       endif
       maxpix = maxval(pixel8)
       minpix = maxval(pixel8)
       if (maxpix > (npixtot-1) .or. minpix < 0) then
          print*,'map constructed from file '//trim(filename)//', with pixels in ',minpix,maxpix
          print*,' wont fit in array with ',npixtot,' elements'
          call fatal_error
       endif
       map(:,:)        = fmiss_effct
       do imap=1,nmaps
          map(pixel8(:),imap) = iqu(:,imap)
       enddo
       imissing = npixtot - obs_npix
       deallocate(iqu)
       deallocate(pixel8)
!     else if (type_fits == 3 .and. (nmaps == 1 .or. nmaps == 3) .and. n_ext == nmaps) then
!        ! cut sky FITS file, reading 1 map (I) or 3 maps (I,Q,U), each from a different extension
    else if (type_fits == 3 .and. (nmaps == 3) .and. n_ext == nmaps) then
       ! cut sky FITS file, reading 3 maps (I,Q,U), each from a different extension
       do imap = 1, nmaps
          extno_f  = extno_i + imap - 1 ! 2016-08-16
          obs_npix = getsize_fits(filename, extno = extno_f)       
          ! one partial map (in binary table with explicit pixel indexing)
          allocate(pixel4(0:obs_npix-1), stat = status)
          allocate(signal(0:obs_npix-1), stat = status)
          call read_fits_cut4(filename, int(obs_npix,kind=i4b), &
               &              pixel4, signal, header=hbuffer, units=units1, extno=extno_f)
          if (present(header) .and. imap == 1) then
             do i=1,nlheader
                header(i) = hbuffer(i)
             enddo
          endif
          if (present(units)) units(imap) = trim(units1)
          maxpix = maxval(pixel4)
          minpix = maxval(pixel4)
          if (maxpix > (npixtot-1) .or. minpix < 0) then
             print*,'map constructed from file '//trim(filename)//', with pixels in ',minpix,maxpix
             print*,' wont fit in array with ',npixtot,' elements'
             call fatal_error
          endif
          map(:,imap)        = fmiss_effct
          map(pixel4(:),imap) = signal(:)
          imissing = npixtot - obs_npix
          deallocate(signal)
          deallocate(pixel4)
       enddo
    else 
       print*,'Unable to read the ',nmaps,' required  map(s) from file '//trim(filename)
       print*,'file type = ',type_fits
       call fatal_error
    endif
    !-----------------------------------------------------------------------
    if (imissing > 0) then
       write(*,'(i12,a,f7.3,a,1pe11.4)') &
            &           imissing,' missing pixels (', &
            &           (100.*imissing)/npixtot,' %),'// &
            &           ' have been set to : ',fmiss_effct
    endif

!    deallocate(unitsm)

    !-----------------------------------------------------------------------
    ! deal with polcconv
    ! print*,'******* in input_map ',nmaps, polar, type_fits, nside_fits

    if ( do_polcconv .and. ( &
         & (nmaps >= 3 .and. type_fits == 3) &  ! cut-sky, multiple maps
         &  .or. &
         & (nmaps >=3 .and. polar==1 .and. type_fits ==2 .and. &
         &    obs_npix==nside2npix(nside_fits)) & ! full-sky, polarized map
         &                  )      ) then 
       if (ipolcconv == 0) then
          print 9000,' The POLCCONV keyword was not found in '//trim(filename)
          print 9000,' COSMO (HEALPix default) will be assumed, and map is unchanged.'
          print 9000,' See HEALPix primer ('//primer_url//') for details.'
       endif
!        if (ipolcconv == 1) then
!           print 9000,' POLCCONV=COSMO found in '//trim(filename)//'. Map is unchanged.'
!        endif
       if (ipolcconv == 2) then
          print 9000,' POLCCONV=IAU found in '//trim(filename)
          map(:,3) = - map(:,3)
          if (present(header)) then
             print 9000,' The sign of the U polarisation is changed in memory,'&
                  & //' and the header is updated.'
             call add_card(header, 'POLCCONV', 'COSMO', &
                  comment='Changed by input_map', update=.true.)
          else
             print 9000,' The sign of the U polarisation is changed in memory.'
          endif
          print 9000,' See HEALPix primer ('//primer_url//') for details.'
       endif
       if (ipolcconv == 3) then
          call get_card(hbuffer,'POLCCONV',polcconv)
          print 9000,' POLCCONV='//trim(polcconv)//' found in '//trim(filename)
          print 9000,' It is neither COSMO nor IAU.   Aborting!'
          print 9000,' See HEALPix primer ('//primer_url//') for details.'
          call fatal_error
       endif
    endif
9000 format(a)
    if (allocated(hbuffer)) deallocate(hbuffer)

    RETURN
  END subroutine input_map8_KLOAD
  !=======================================================================
  !     Read a FITS file
  !     This routine is used for reading MAPS by anafast.
  !     modified Feb 03 for units argument to run with Sun compiler
  !     Jan 2005, EH, improved for faster writing
  !=======================================================================
#ifndef NO64BITS
  subroutine read_bintab4_KLOAD(filename, map, npixtot, nmaps, nullval, anynull, header, units, extno)
    character(len=*),                          intent(IN)  :: filename
    integer(I4B),                              intent(IN)  :: npixtot
    integer(I4B),                              intent(IN)  :: nmaps
    real(KMAP),      dimension(0:,1:),         intent(OUT) :: map
    real(KMAP),                                intent(OUT) :: nullval
    logical(LGT),                              intent(OUT) :: anynull
    character(LEN=*), dimension(1:), optional, intent(OUT) :: header
    character(LEN=*), dimension(1:), optional, intent(OUT) :: units
    integer(I4B)                   , optional, intent(IN) :: extno

    integer(i8b)::  npixtot8

    npixtot8 = int(npixtot,kind=i8b)
    call read_bintab8_KLOAD(filename, map, npixtot8, nmaps, nullval, anynull, header, units, extno)
    return

  end subroutine read_bintab4_KLOAD
#endif
  subroutine read_bintab8_KLOAD(filename, map, npixtot, nmaps, nullval, anynull, header, units, extno)
  !=======================================================================
    character(len=*),                          intent(IN)  :: filename
    integer(I8B),                              intent(IN)  :: npixtot
    integer(I4B),                              intent(IN)  :: nmaps
    real(KMAP),      dimension(0:,1:),         intent(OUT) :: map
    real(KMAP),                                intent(OUT) :: nullval
    logical(LGT),                              intent(OUT) :: anynull
    character(LEN=*), dimension(1:), optional, intent(OUT) :: header
    character(LEN=*), dimension(1:), optional, intent(OUT) :: units
    integer(I4B)                   , optional, intent(IN) :: extno

    integer(I4B) :: nl_header, len_header, nl_units, len_units
    integer(I4B) :: status,unit,readwrite,blocksize,naxes(2),nfound, naxis
    integer(I4B) :: group, firstpix, i, npix32
    real(KMAP)   :: blank, testval
    real(DP)     :: bscale,bzero
    character(len=80) :: comment
    logical(LGT) :: extend
    integer(I4B) :: nmove, hdutype, hdunum
    integer(I4B) :: frow, imap
    integer(I4B) :: datacode, width
    LOGICAL(LGT) ::  anynull_i

    integer(I4B),     parameter            :: MAXDIM = MAXDIM_TOP !number of columns in the extension
    integer(i8b)                      :: npix_old
    integer(i8b), dimension(1:MAXDIM) :: npix
    integer(i8b), dimension(1:MAXDIM) :: i0, i1
    integer(i4b), dimension(1:MAXDIM) :: repeat
    integer(i4b)                      :: nrow2read, nelem

    integer(I4B)                           :: nrows, tfields, varidat
    character(len=20), dimension(1:MAXDIM) :: ttype, tform, tunit
    character(len=20)                      :: extname
    character(len=*), parameter            :: code='read_bintab'
    !-----------------------------------------------------------------------
    status=0

    unit = 146
    naxes(1) = 1
    naxes(2) = 1
    nfound = -1
    anynull = .false.
    bscale = 1.0d0
    bzero = 0.0d0
    blank = -2.e25
    nullval = bscale*blank + bzero
    comment=''
    ttype=''
    tform=''
    tunit=''
    extname=''

    nl_header = 0
    if (present(header)) then
       nl_header = size(header)
       len_header = 80
    endif

    nl_units = 0
    if (present(units)) then
       nl_units = size(units)
       len_units = min(80,len(units(1))) ! due to SUN compiler bug
    endif

    readwrite=0
    !call ftopen(unit,filename,readwrite,blocksize,status)
    call ftnopn(unit,filename,readwrite, status) !open primary HDU or specified HDU
    if (status > 0) call printerror(status)
    !     -----------------------------------------
    call ftghdn(unit, hdunum)

    if (hdunum == 1) then  ! in primary HDU
       !     determines the presence of image
       call ftgkyj(unit,'NAXIS', naxis, comment, status)
       if (status > 0) call printerror(status)

       !     determines the presence of an extension
       call ftgkyl(unit,'EXTEND', extend, comment, status)
       if (status > 0) then
          extend = .false.
          status = 0 ! no extension : 
       !     to be compatible with first version of the code
       endif
    endif

    if (naxis > 0 .and. .not.extend .and. hdunum==1) then ! there is an image
       !        determine the size of the image (look naxis1 and naxis2)
       call ftgknj(unit,'NAXIS',1_i4b,2_i4b,naxes,nfound,status)

       !        check that it found only NAXIS1
       if (nfound == 2 .and. naxes(2) > 1) then
          print *,'multi-dimensional image'
          print *,'expected 1-D data.'
          call fatal_error
       end if

       if (nfound < 1) then
          call printerror(status)
          print *,'can not find NAXIS1.'
          call fatal_error
       endif

       npix(1)=naxes(1)
       if (npix(1) /= npixtot) then
          print *,'WARNING: found ',npix(1),' pixels in '//trim(filename)
          print *,'         expected ',npixtot
          npix(1) = min(npix(1), npixtot)
          print *,'         only ',npix(1),' will be read'
       endif

       call ftgkyd(unit,'BSCALE',bscale,comment,status)
       if (status == 202) then ! BSCALE not found
          bscale = 1.0d0
          status = 0
       endif
       call ftgkyd(unit,'BZERO', bzero, comment,status)
       if (status == 202) then ! BZERO not found
          bzero = 0.0d0
          status = 0
       endif
       call f90ftgky_(unit, 'BLANK', blank, comment, status)
       if (status == 202) then ! BLANK not found 
          ! (according to fitsio BLANK is integer)
          blank = -2.e25
          status = 0
       endif
       nullval = bscale*blank + bzero

       !        -----------------------------------------

       group=1
       firstpix = 1
       npix32 = npix(1)
       call f90ftgpv_(unit, group, firstpix, npix32, nullval, map(0:npix(1)-1,1), anynull, status)
       ! if there are any NaN pixels, (real data)
       ! or BLANK pixels (integer data) they will take nullval value
       ! and anynull will switch to .true.
       ! otherwise, switch it by hand if necessary
       testval = 1.e-6 * ABS(nullval)
       do i=0, npix(1)-1
          if (ABS(map(i,1)-nullval) < testval) then
             anynull = .true.
             goto 111
          endif
       enddo
111    continue

    else if (extend .or. hdunum>1) then
       if (hdunum == 1) then
          nmove = +1
          if (present(extno)) nmove = +1 + extno
          call ftmrhd(unit, nmove, hdutype, status)
       else
          call ftghdt(unit, hdutype, status)
       endif

       call assert (hdutype==2, 'this is not a binary table')

       !        reads all the keywords
       call ftghbn(unit, MAXDIM, &
            &        nrows, tfields, ttype, tform, tunit, extname, varidat, &
            &        status)

       if (tfields < nmaps) then
          print *,'found ',tfields,' maps in file '//trim(filename)
          print *,'expected ',nmaps
          call fatal_error
       endif

       !        finds the bad data value
       call f90ftgky_(unit, 'BAD_DATA', nullval, comment, status)
       if (status == 202) then ! bad_data not found
          if (KMAP == SP) nullval = s_bad_value ! default value
          if (KMAP == DP) nullval = d_bad_value ! default value
          status = 0
       endif

       if (nl_header > 0) then
          do i=1,nl_header
             header(i)(1:len_header) = ""
          enddo
          call get_clean_header(unit, header, filename, status)
          status = 0
       endif

       if (nl_units > 0) then
          do i=1,nl_units
             units(i)(1:len_units) = 'unknown' ! default
          enddo
          do imap = 1, min(nmaps, nl_units)
             units(imap)(1:len_units) = adjustl(tunit(imap))
          enddo
       endif

       npix_old = npixtot
       do imap = 1, nmaps
          !parse TFORM keyword to find out the length of the column vector
          call ftbnfm(tform(imap), datacode, repeat(imap), width, status)
          npix(imap) = int(nrows,i8b) * repeat(imap)
          if (npix(imap) /= npixtot .and. npix_old /= npix(imap)) then
             print *,'WARNING: found ',npix(imap),' pixels in '//trim(filename)//', column ',imap
             print *,'         expected ',npixtot,' or ',npix_old
             npix_old = npix(imap)
             npix(imap) = min(npix(imap), npixtot)
             print *,'         only  ',npix(imap),' will be read'
          endif
       enddo


       call ftgrsz(unit, nrow2read, status)
       nrow2read = max(nrow2read, 1)
       firstpix  = 1  ! starting position in FITS within row, 1 based
       i0(:) = 0_i8b  ! starting element in array, 0 based
       do frow = 1, nrows, nrow2read
          do imap = 1, nmaps
             i1(imap) = min(i0(imap) + int(nrow2read,i8b) * repeat(imap), npix(imap)) - 1_i8b
             nelem = i1(imap) - i0(imap) + 1
             call f90ftgcv_(unit, imap, frow, firstpix, nelem, &
                  & nullval, map(i0(imap):i1(imap),imap), anynull_i, status)
             anynull = anynull .or. anynull_i
             i0(imap) = i1(imap) + 1_i8b
          enddo
       enddo
       ! sanity check
       do imap = 1, nmaps
          if (i0(imap) /= npix(imap)) then
             call fatal_error('something wrong during piece wise reading')
          endif
       enddo

    else ! no image no extension, you are dead, man
       call fatal_error(' No image, no extension in '//trim(filename))
    endif
    !     close the file
    call ftclos(unit, status)

    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)

    return
  end subroutine read_bintab8_KLOAD

  !=======================================================================
  subroutine read_conbintab_KLOAD(filename, alms, nalms, units, extno)
    !=======================================================================
    !     Read a FITS file containing alms values
    !
    !     slightly modified to deal with vector column 
    !     in binary table       EH/IAP/Jan-98
    !
    !     Used by synfast when reading a binary file with alms for cons.real.
    !                        FKH/Apr-99
    !
    !   extno : optional, number of extension to be read, default=0: first extension
    !     Jan 2005, EH, improved for faster reading
    !=======================================================================
    CHARACTER(LEN=*),                   INTENT(IN) :: filename
    INTEGER(I4B),                       INTENT(IN) :: nalms !, nlheader
    REAL(KMAP), DIMENSION(0:nalms-1,1:6), INTENT(OUT) :: alms
    CHARACTER(LEN=*), dimension(1:), optional,  INTENT(OUT) :: units
    INTEGER(I4B)                   , optional,  INTENT(IN) :: extno

    REAL(KMAP)                                        :: nullval
    LOGICAL(LGT)                                    ::  anynull

    INTEGER(I4B), DIMENSION(:), allocatable :: lm
    INTEGER(I4B) :: status,unit,readwrite,blocksize,naxes(2),nfound, naxis
    INTEGER(I4B) :: npix
    CHARACTER(LEN=80) :: comment
    LOGICAL(LGT) :: extend
    INTEGER(I4B) :: nmove, hdutype
    INTEGER(I4B) :: frow, imap
    INTEGER(I4B) :: datacode, repeat, width
    integer(I4B) :: i, l, m
    integer(i4b) :: nrow2read, nelem
    integer(i8b) :: i0, i1

    INTEGER(I4B), PARAMETER :: MAXDIM = MAXDIM_TOP !number of columns in the extension
    INTEGER(I4B) :: nrows, tfields, varidat
    CHARACTER(LEN=20), dimension(1:MAXDIM) :: ttype, tform, tunit
    CHARACTER(LEN=20) :: extname
    character(len=*), parameter :: code="read_conbintab"

    !-----------------------------------------------------------------------
    status=0

    unit = 145
    naxes(1) = 1
    naxes(2) = 1
    nfound = -1
    anynull = .false.
    alms=0.  ! set output to 0.
    readwrite=0
    comment=''
    ttype=''
    tform=''
    tunit=''
    extname=''
    call ftopen(unit,filename,readwrite,blocksize,status)
    if (status > 0) call printerror(status)
    !     -----------------------------------------

    !     determines the presence of image
    call ftgkyj(unit,'NAXIS', naxis, comment, status)
    if (status > 0) call printerror(status)

    !     determines the presence of an extension
    call ftgkyl(unit,'EXTEND', extend, comment, status)
    if (status > 0) status = 0 ! no extension : 
    !     to be compatible with first version of the code

    if (.not.extend) then
       print*,'No extension!'
       call fatal_error
    endif

    ! go to assigned extension (if it exists)
    nmove = +1
    if (present(extno)) nmove = +1 + extno
    call ftmrhd(unit, nmove, hdutype, status)
    if (status > 0) then
       ! if required extension not found: 
       ! print a warning, fill with dummy values, return to calling routine
       print*,code//' WARNING: the extension ',extno,' was not found in ',trim(filename)
       alms(0:nalms-1,1)=0.    ! l = 0
       alms(0:nalms-1,2)=1.    ! m = 1
       alms(0:nalms-1,3:6)=0.
       status = 0
       call ftclos(unit, status)        !     close the file
       return
    endif

    if (hdutype /= 2) then ! not a binary table
       print*, 'this is not a binary table'
       call fatal_error
    endif

    !        reads all the keywords
    call ftghbn(unit, MAXDIM, &
         &        nrows, tfields, ttype, tform, tunit, extname, varidat, &
         &        status)

    if ((tfields/=3).and.(tfields/=5)) then
       print *,'found ',tfields,' columns in the file'
       print *,'expected 3 or 5'
       call fatal_error
    endif
    !        finds the bad data value
!     if (KMAP == SP) call ftgkye(unit,'BAD_DATA',nullval,comment,status)
!     if (KMAP == DP) call ftgkyd(unit,'BAD_DATA',nullval,comment,status)
    call f90ftgky_(unit, 'BAD_DATA', nullval, comment, status)
    if (status == 202) then ! bad_data not found
       if (KMAP == SP) nullval = s_bad_value ! default value
       if (KMAP == DP) nullval = d_bad_value ! default value
       status = 0
    endif

    !          if (nlheader > 0) then
    !             header = ""
    !             status = 0
    !             call get_clean_header(unit, header, filename, status)
    !          endif

    if (present(units)) then
       units(1) = tunit(2)
    endif

    !parse TFORM keyword to find out the length of the column vector
    call ftbnfm(tform(1), datacode, repeat, width, status)
    npix = nrows * repeat
    if (npix /= nalms) then
       print *,'found ',npix,' alms'
       print *,'expected ',nalms
       call fatal_error
    endif

    call ftgrsz(unit, nrow2read, status)
    nrow2read = max(nrow2read, 1)
    nelem = nrow2read * repeat
    i0 = 0_i8b
    allocate(lm(0:nelem-1))
    do frow = 1, nrows, nrow2read
       i1 = min(i0 + nrow2read * repeat, int(npix,i8b)) - 1_i8b
       nelem = i1 - i0 + 1
       ! first column -> index
       call ftgcvj(unit, 1_i4b, frow, 1_i4b, nelem, i_bad_value, &
         &        lm(0), anynull, status)
       call assert(.not. anynull, 'There are undefined values in the table!')
       ! other columns -> a(l,m)
       do imap = 2, tfields
          call f90ftgcv_(unit, imap, frow, 1_i4b, nelem, nullval, &
               &        alms(i0:i1,imap+1), anynull, status)
          call assert (.not. anynull, 'There are undefined values in the table!')
       enddo
       ! recoding of the mapping, EH, 2002-08-08
       ! lm = l*l + l + m + 1
       do i = i0, i1
          l = int(sqrt(   real(lm(i-i0)-1, kind = DP)  ) )
          m = lm(i-i0) - l*(l+1) - 1
          ! check that round-off did not screw up mapping
          if (abs(m) > l) then
             print*,'Inconsistent l^2+l+m+1 -> l,m mapping'
             print*,l, m, l*(l+1)+m+1, lm(i-i0)
             call fatal_error
          endif
          alms(i,1) = real( l, kind=KMAP)
          alms(i,2) = real( m, kind=KMAP)
       enddo
       i0 = i1 + 1_i8b
    enddo
    deallocate(lm)
    ! sanity check
    if (i0 /= npix) then
       print*,'something wrong during piece-wise reading'
       call fatal_error
    endif

!     !reads the columns
!     ! first column : index -> (l,m)
!     allocate(lm(0:nalms-1))
!     column = 1
!     frow = 1
!     firstpix = 1
!     npix = nrows * repeat
!     if (npix /= nalms) then
!        print *,'found ',npix,' alms'
!        print *,'expected ',nalms
!        call fatal_error
!     endif
!     call ftgcvj(unit, column, frow, firstpix, npix, i_bad_value, &
!          &        lm(0), anynull, status)
!     call assert (.not. anynull, 'There are undefined values in the table!')
!     ! recoding of the mapping, EH, 2002-08-08
!     ! lm = l*l + l + m + 1
!     do i = 0, nalms - 1
!        l = int(sqrt(   real(lm(i)-1, kind = DP)  ) )
!        m = lm(i) - l*(l+1) - 1
!        ! check that round-off did not screw up mapping
!        if (abs(m) > l) then
!           print*,'Inconsistent l^2+l+m+1 -> l,m mapping'
!           print*,l, m, l*(l+1)+m+1, lm(i)
!           call fatal_error
!        endif
!        alms(i,1) = real( l, kind=KMAP )
!        alms(i,2) = real( m, kind=KMAP )
!     enddo
!     deallocate(lm)
!
!     do imap = 2, tfields
!        !parse TFORM keyword to find out the length of the column vector
!        call ftbnfm(tform(imap), datacode, repeat, width, status)
!
!        !reads the columns
!        column = imap
!        frow = 1
!        firstpix = 1
!        npix = nrows * repeat
!        if (npix /= nalms) then
!           print *,'found ',npix,' alms'
!           print *,'expected ',nalms
!           call fatal_error
!        endif
!        call f90ftgcv_(unit, column, frow, firstpix, npix, nullval, &
!             &        alms(0:npix-1,imap+1), anynull, status)
!        call assert (.not. anynull, 'There are undefined values in the table!')
!     enddo
    !     close the file
    call ftclos(unit, status)


    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)
    return
  end subroutine read_conbintab_KLOAD

#ifndef NO64BITS
  subroutine write_bintab4_KLOAD(map, npix, nmap, header, nlheader, filename, extno)
  !=======================================================================
    INTEGER(I4B),     INTENT(IN) :: npix, nmap, nlheader
    REAL(KMAP),       INTENT(IN), DIMENSION(0:npix-1,1:nmap) :: map
    CHARACTER(LEN=*), INTENT(IN), DIMENSION(1:nlheader)      :: header
    CHARACTER(LEN=*), INTENT(IN)           :: filename
    INTEGER(I4B)    , INTENT(IN), optional :: extno

    integer(i8b) :: npix8

    npix8 = npix
    if (present(extno)) then
       call write_bintab8_KLOAD(map, npix8, nmap, header, nlheader, filename, extno)
    else
       call write_bintab8_KLOAD(map, npix8, nmap, header, nlheader, filename)
    endif
    return
  end subroutine write_bintab4_KLOAD
#endif

  subroutine write_bintab8_KLOAD(map, npix, nmap, header, nlheader, filename, extno)
    !=======================================================================
    !     Create a FITS file containing a binary table extension with 
    !     the temperature map in the first column
    !     written by EH from writeimage and writebintable 
    !     (fitsio cookbook package)
    !
    !     slightly modified to deal with vector column (ie TFORMi = '1024E')
    !     in binary table       EH/IAP/Jan-98
    !
    !     simplified the calling sequence, the header sould be filled in
    !     before calling the routine
    !
    !     July 21, 2004: SP version
    !     Jan 2005, EH, improved for faster writing
    !     2009-02-25: EH, accepts I4B and I8B npix
    !=======================================================================
    use pix_tools, only: npix2nside
    INTEGER(I8B),     INTENT(IN) :: npix
    INTEGER(I4B),     INTENT(IN) :: nmap, nlheader
!     REAL(KMAP),       INTENT(IN), DIMENSION(0:npix-1,1:nmap), target :: map
    REAL(KMAP),       INTENT(IN), DIMENSION(0:npix-1,1:nmap) :: map
    CHARACTER(LEN=*), INTENT(IN), DIMENSION(1:nlheader)      :: header
    CHARACTER(LEN=*), INTENT(IN)           :: filename
    INTEGER(I4B)    , INTENT(IN), optional :: extno

    INTEGER(I4B) ::  status,unit,blocksize,bitpix,naxis,naxes(1)
    INTEGER(I4B) ::  i, nside
    LOGICAL(LGT) ::  simple,extend
    CHARACTER(LEN=80) :: comment, ch

    INTEGER(I4B), PARAMETER :: MAXDIM = MAXDIM_TOP !number of columns in the extension
    INTEGER(I4B) :: nrows, tfields, varidat
    integer(i4b) :: repeat, nrow2write, nelem
    integer(i8b) :: i0, i1
    INTEGER(I4B) :: frow,  felem, colnum, hdutype
    CHARACTER(LEN=20) :: ttype(MAXDIM), tform(MAXDIM), tunit(MAXDIM), extname
    CHARACTER(LEN=10) ::  card
    CHARACTER(LEN=2) :: stn
    INTEGER(I4B) :: itn, extno_i
    character(len=filenamelen) sfilename
    character(len=1)          :: pform 
!    character(len=80)         :: junk_k, junk_v, junk_c
    !-----------------------------------------------------------------------

    if (KMAP == SP) pform = 'E'
    if (KMAP == DP) pform = 'D'
    status=0
    unit = 144
    blocksize=1

    extno_i = 0
    if (present(extno)) extno_i = extno

    if (extno_i == 0) then
       !*************************************
       !     create the new empty FITS file
       !*************************************
       call ftinit(unit,filename,blocksize,status)
       if (status > 0) call fatal_error("Error while creating file " &
            & //trim(filename) &
            & //". Check path and/or access rights.")

       !     -----------------------------------------------------
       !     initialize parameters about the FITS image
       simple=.true.
       bitpix=32     ! integer*4
       naxis=0       ! no image
       naxes(1)=0
       extend=.true. ! there is an extension

       !     ----------------------
       !     primary header
       !     ----------------------
       !     write the required header keywords
       call ftphpr(unit,simple,bitpix,naxis,naxes,0_i4b,1_i4b,extend,status)

       !     writes supplementary keywords : none

       !     write the current date
       call ftpdat(unit,status) ! format (yyyy-mm-dd)

       !     ----------------------
       !     image : none
       !     ----------------------

       !     ----------------------
       !     extension
       !     ----------------------
    else

       !*********************************************
       !     reopen an existing file and go to the end
       !*********************************************
       ! remove the leading '!' (if any) when reopening the same file
       sfilename = adjustl(filename)
       if (sfilename(1:1) == '!') sfilename = sfilename(2:filenamelen)
       call ftopen(unit,sfilename,1_i4b,blocksize,status)
       call ftmahd(unit,1_i4b+extno_i,hdutype,status)

    endif

    !     creates an extension
    call ftcrhd(unit, status)

    !     writes required keywords
    tfields  = nmap
    repeat   = 1024
    if (npix < 1024) repeat = 1 ! for nside <= 8
    ! for large npix increase repeat so that nrows < 2^31-1
    nside = npix2nside(npix)
    if (nside > 1024*256) repeat = nside/256
    nrows    = npix / repeat ! naxis1
    ch = string(repeat, format='(i8)')
    ch = trim(adjustl(ch))//pform
    tform(1:nmap) = ch
    ttype(1:nmap) = 'simulation'   ! will be updated
    tunit(1:nmap) = ''      ! optional, will not appear
    extname  = ''      ! optional, will not appear
    varidat  = 0
    call ftphbn(unit, nrows, tfields, ttype, tform, tunit, &
         &     extname, varidat, status)

    !     write the header literally, putting TFORM1 and TUNIT1 at the desired place
    if (KMAP == SP) comment = 'data format of field: 4-byte REAL'
    if (KMAP == DP) comment = 'data format of field: 8-byte REAL'
    do i=1,nlheader
       card = header(i)
       if (card(1:5) == 'TTYPE') then ! if TTYPEi is explicitely given
          stn = card(6:8)
          read(stn,'(i3)') itn
          ! discard at their original location:
          call ftdkey(unit,'TTYPE'//stn,status)  ! old TTYPEi and
          status = 0
          call ftdkey(unit,'TFORM'//stn,status)  !     TFORMi
          status = 0
          if (itn <= tfields) then ! only put relevant information 2008-08-27
             call putrec(unit,header(i), status)    ! write new TTYPE1
             status = 0
             call ftpkys(unit,'TFORM'//stn,tform(itn),comment,status) ! and write new TFORM1 right after
          endif
       elseif (header(i)/=' ') then
          call putrec(unit,header(i), status)
       endif
       status = 0
    enddo

    !     write the extension buffer by buffer
    call ftgrsz(unit, nrow2write, status)
    nrow2write = max(nrow2write, 1)
    felem  = 1  ! starting position in FITS (element), 1 based
    i0 = 0_i8b  ! starting element in array, 0 based
    do frow = 1, nrows, nrow2write
       i1 = min(i0 + nrow2write * repeat, npix) - 1_i8b
       nelem = i1 - i0 + 1
       do colnum = 1, nmap
          call f90ftpcl_(unit, colnum, frow, felem, nelem, map(i0:i1, colnum), status)
       enddo
       i0 = i1 + 1_i8b
    enddo
    ! sanity check
    if (i0 /= npix) then
       call fatal_error("something wrong during piece wise writing")
    endif

    !     ----------------------
    !     close and exit
    !     ----------------------

    !     close the file and free the unit number
    call ftclos(unit, status)

    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)

    return
  end subroutine write_bintab8_KLOAD
  !=======================================================================
  subroutine write_asctab_KLOAD &
       &  (clout, lmax, ncl, header, nlheader, filename, extno)
    !=======================================================================
    !     Create a FITS file containing an ASCII table extension with 
    !     the measured power spectra
    !     written by EH from writeimage and writeascii 
    !     (fitsio cookbook package)
    !
    !     clout = power spectra with l in (0:lmax)
    !     ncl = number of spectra
    !     header = FITS header to be put on top of the file
    !     nlheader = number of lines of the header
    !     filename = FITS output file name
    !
    !     2015-09-10: EH, added extno
    !=======================================================================
    !
    INTEGER(I4B),      INTENT(IN)           :: lmax, ncl,nlheader
    REAL(KMAP),        INTENT(IN)           ::  clout(0:lmax,1:ncl)
    CHARACTER(LEN=80), INTENT(IN), DIMENSION(1:nlheader) :: header
    CHARACTER(LEN=*),  INTENT(IN)           :: filename
    INTEGER(I4B)    ,  INTENT(IN), optional :: extno

    INTEGER(I4B) :: bitpix,naxis,naxes(1)
    LOGICAL(LGT) :: simple,extend
    CHARACTER(LEN=10) ::  card

    INTEGER(I4B), PARAMETER :: MAXDIM = MAXDIM_TOP
    INTEGER(I4B) ::  status,unit,blocksize,tfields,nrows,rowlen
    INTEGER(I4B) ::  nspace,tbcol(MAXDIM),colnum,frow,felem,hdutype
    CHARACTER(LEN=16) :: ttype(MAXDIM),tform(MAXDIM),tunit(MAXDIM),extname
    CHARACTER(LEN=80) :: comment, card_tbcol
    CHARACTER(LEN=2) :: stn
    INTEGER(I4B) :: itn, i, extno_i
    character(len=filenamelen) sfilename
    character(len=6)        :: form
    !=======================================================================
    if (KMAP == SP) form = 'E15.7'
    if (KMAP == DP) form = 'D24.15'
    status=0
    unit = 109
    blocksize=1

    extno_i = 0
    if (present(extno)) extno_i = extno

    if (extno_i == 0) then
       call ftinit(unit,filename,blocksize,status)
       if (status > 0) call fatal_error("Error while creating file " &
            & //trim(filename) &
            & //". Check path and/or access rights.")
       
       !     -----------------------------------------------------
       !     initialize parameters about the FITS image
       simple=.true.
       bitpix=32     ! integer*4
       naxis=0       ! no image
       naxes(1)=0
       extend=.true. ! there is an extension
       
       !     ----------------------
       !     primary header
       !     ----------------------
       !     write the required header keywords
       call ftphpr(unit,simple,bitpix,naxis,naxes,0_i4b,1_i4b,extend,status)
       
       !     writes supplementary keywords : none
       
       !     write the current date
       call ftpdat(unit,status) ! (format ccyy-mm-dd)
       
       !     ----------------------
       !     image : none
       !     ----------------------
       
       !     ----------------------
       !     extension
       !     ----------------------
    else
       !*********************************************
       !     reopen an existing file and go to the end
       !*********************************************
       ! remove the leading '!' (if any) when reopening the same file
       sfilename = adjustl(filename)
       if (sfilename(1:1) == '!') sfilename = sfilename(2:filenamelen)
       call ftopen(unit,sfilename,1_i4b,blocksize,status)
       call ftmahd(unit,1_i4b+extno_i,hdutype,status)

    endif

    !     append a new empty extension onto the end of the primary array
    call ftcrhd(unit,status)

    !     define parameters for the ASCII table
    nrows   = lmax+1
    tfields = ncl
    tform(1:ncl) = form
    ttype(1:ncl) = 'power spectrum' ! is updated by the value given in the header
    tunit(1:ncl) = '' ! optional, will not appear
    extname      = '' ! optional, will not appear

    !     calculate the starting position of each column, and the total row length
    nspace=1
    call ftgabc(tfields,tform,nspace,rowlen,tbcol,status)

    !     write the required header parameters for the ASCII table
    call ftphtb(unit,rowlen,nrows,tfields,ttype,tbcol,tform,tunit, &
         &            extname,status)

    !     write the header literally, putting TFORM1 at the desired place
    comment = ''
    card_tbcol=''
    do i=1,nlheader
       card = header(i)
       if (card(1:5) == 'TTYPE') then ! if TTYPE1 is explicitely given
          stn = card(6:8)
          read(stn,'(i3)') itn
          ! discard at their original location:
          !!!!!!!call ftdkey(unit,card(1:6),status)  !         old TTYPEi 2019-10-14
          call ftdkey(unit,'TTYPE'//stn,status)  !         old TTYPEi
          status = 0
          call ftdkey(unit,'TFORM'//stn,status) !           TFORMi
          status = 0
          call ftgcrd(unit,'TBCOL'//stn,card_tbcol,status)
          status = 0
          call ftdkey(unit,'TBCOL'//stn,status) !           TBCOLi
          status = 0
          ! and rewrite
          if (itn <= tfields) then
             call putrec(unit,card_tbcol,status) !             TBCOLi
             status = 0
             call ftpkys(unit,'TFORM'//stn,tform(itn),comment,status) ! TFORMi right after
             status = 0
             call putrec(unit,header(i), status)   !           TTYPEi
          endif
       elseif (header(i)/=' ') then
          call putrec(unit,header(i), status)
       endif
       status = 0
    enddo

    frow=1
    felem=1
    do colnum = 1, ncl
       call f90ftpcl_(unit, colnum, frow, felem, nrows, clout(0:nrows-1,colnum), status)  
    enddo

    !     close the FITS file
    call ftclos(unit, status)

    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)

    return
  end subroutine write_asctab_KLOAD
  !=======================================================================
  !   DUMP_ALMS
  !
  !     Create/extend a FITS file containing a binary table extension with 
  !     the a_lm coefficients in each extension.
  !
  !     Format of alm FITS file: in each extension, 3 columns: 
  !            index=l*l+l+m+1, real(a_lm), imag(a_lm)
  !      the a_lm are obtained using the complex Y_lm
  !      only the modes with m>=0 are stored (because the map is assumed real)
  !
  !     The extensions contain, in general, T, E (= G) and B (= C) in that order
  !     First extension has extno = 0
  !
  !     Adapted from write_bintab, FKH/Apr-99
  !   SP/DP overloaded, 2004-12, EH
  !    reduced size of intermediate arrays (alms_out, lm)
  !
  subroutine dump_alms_KLOAD(filename, alms, nlmax, header, nlheader, extno)
    !=======================================================================
    INTEGER(I4B),      INTENT(IN) :: nlmax, nlheader, extno
    COMPLEX(KMAPC),      INTENT(IN), DIMENSION(0:,0:) :: alms
    CHARACTER(LEN=80), INTENT(IN), DIMENSION(1:nlheader) :: header
    CHARACTER(LEN=*),  INTENT(IN)               :: filename
    ! local variables
    INTEGER(I4B),   DIMENSION(:),   allocatable :: lm
    REAL(KMAP),     DIMENSION(:,:), allocatable :: alms_out
    INTEGER(I4B) ::  status,unit,blocksize,bitpix,naxis,naxes(1)
    INTEGER(I4B) ::  i,l,m,cnt,hdutype, nmmax
    LOGICAL(LGT) ::  simple,extend, found_lmax, found_mmax
    CHARACTER(LEN=80) :: comment

    integer(i8b) :: npix
    INTEGER(I4B), PARAMETER :: MAXDIM = MAXDIM_TOP !number of columns in the extension
    INTEGER(I4B) :: nrows, tfields, varidat
    INTEGER(I4B) :: frow,  felem
    CHARACTER(LEN=20) :: ttype(MAXDIM), tform(MAXDIM), tunit(MAXDIM), extname
    CHARACTER(LEN=10) ::  card
    CHARACTER(LEN=2) :: stn
    INTEGER(I4B) :: itn
    character(len=filenamelen) sfilename
    character(len=1) :: pform
    character(len=*), parameter :: code = 'dump_alms'
    !-----------------------------------------------------------------------

    if (KMAP == SP) pform = 'E'
    if (KMAP == DP) pform = 'D'
    
    nmmax = size(alms,2) - 1
    if (nmmax < 0 .or. nmmax > nlmax) call fatal_error('inconsistent Lmax and Mmax in dump_alms')
!!    npix=((nlmax+1)*(nlmax+2))/2
    npix = ((nmmax+1_i8b)*(2_i8b*nlmax-nmmax+2))/2
    if (npix > MAX_I4B .or. nlmax >= int(sqrt(real(MAX_I4B,kind=dp)))) then
       print*,code//'> Index of a_lm too large for current file format,'
       print*,code//'> '//trim(filename)//' not written.'
       call fatal_error
    endif

    status=0
    unit = 141
    blocksize=1

    if (extno==0) then
       !*********************************************
       !     create the new empty FITS file
       !*********************************************
       call ftinit(unit,filename,blocksize,status)
       if (status > 0) call fatal_error("Error while creating file " &
            & //trim(filename) &
            & //". Check path and/or access rights.")

       !     -----------------------------------------------------
       !     initialize parameters about the FITS image
       simple=.true.
       bitpix=32     ! integer*4
       naxis=0       ! no image
       naxes(1)=0
       extend=.true. ! there is an extension

       !     ----------------------
       !     primary header
       !     ----------------------
       !     write the required header keywords
       call ftphpr(unit,simple,bitpix,naxis,naxes,0_i4b,1_i4b,extend,status)

       !     writes supplementary keywords : none

       !     write the current date
       call ftpdat(unit,status) ! format (yyyy-mm-dd)

       !     ----------------------
       !     image : none
       !     extension
       !     ----------------------
    else

       !*********************************************
       !     reopen an existing file and go to the end
       !*********************************************
       ! remove the leading '!' (if any) when reopening the same file
       sfilename = adjustl(filename)
       if (sfilename(1:1) == '!') sfilename = sfilename(2:filenamelen)
       call ftopen(unit,sfilename,1_i4b,blocksize,status)
       call ftmahd(unit,1_i4b+extno,hdutype,status)

    endif

    !     creates an extension
    call ftcrhd(unit, status)

    !     writes required keywords
    nrows    = npix  ! naxis1
    tfields  = 3
    tform = '' ; ttype = '' ; tunit = ''
    tform(1)(1:2)=  '1J' ! necessary with ifort 11.1.076
    tform(2)(1:2) = '1'//pform
    tform(3)(1:2) = '1'//pform
    ttype(1)(1:15) = 'index=l^2+l+m+1'
    ttype(2)(1:15) = 'alm (real)     '
    ttype(3)(1:15) = 'alm (imaginary)' 
    tunit(1:3) = ''      ! optional, will not appear
    extname  = ''      ! optional, will not appear
    varidat  = 0
    call ftphbn(unit, nrows, tfields, ttype, tform, tunit, &
         &     extname, varidat, status)

    !     write the header literally, putting TFORM1 at the desired place
    do i=1,nlheader
       card = header(i)
       if (card(1:5) == 'TTYPE') then ! if TTYPEi is explicitely given
          stn = card(6:8)
          read(stn,'(i3)') itn
          ! discard at their original location:
          call ftdkey(unit,'TTYPE'//stn,status)  ! old TTYPEi and  ! remove
          status = 0
          call ftdkey(unit,'TFORM'//stn,status)  !     TFORMi
          status = 0
          if (itn <= tfields) then ! only put relevant information 2008-08-27
             call putrec(unit,header(i), status)           ! write new TTYPEi
             status = 0
             if (itn==1) then
                comment = 'data format of field: 4-byte INTEGER'
             else
                if (KMAP == SP) comment = 'data format of field: 4-byte REAL'
                if (KMAP == DP) comment = 'data format of field: 8-byte REAL'
             endif
             call ftpkys(unit,'TFORM'//stn,tform(itn),comment,status) ! and write new TFORM1 right after
          endif
       elseif (header(i)/=' ') then
          call putrec(unit,header(i), status)
       endif
       status = 0
    enddo

    call ftukyj(unit, 'MIN-LPOL', 0_i4b, 'Minimum L multipole order',  status)
    call ftukyj(unit, 'MAX-LPOL', nlmax, 'Maximum L multipole order',  status)
    call ftukyj(unit, 'MAX-MPOL', nmmax, 'Maximum M multipole degree', status)

    allocate(lm      (0:nmmax))
    allocate(alms_out(0:nmmax,1:2))
    !     write the extension one column by one column
    frow   = 1  ! starting position (row)
    felem  = 1  ! starting position (element)
    do l = 0, nlmax
       cnt = 0
       do m = 0, min(l,nmmax)
          lm(cnt) = l**2 + l + m + 1
          alms_out(cnt,1)=REAL( alms(l,m))
          alms_out(cnt,2)=AIMAG(alms(l,m))
          cnt = cnt + 1
       enddo
       call ftpclj(   unit, 1_i4b, frow, felem, cnt, lm,                  status)
       call f90ftpcl_(unit, 2_i4b, frow, felem, cnt, alms_out(0:cnt-1,1), status)
       call f90ftpcl_(unit, 3_i4b, frow, felem, cnt, alms_out(0:cnt-1,2), status)
       frow = frow + cnt
    enddo
    deallocate(lm)
    deallocate(alms_out)

    !     ----------------------
    !     close and exit
    !     ----------------------

    !     close the file and free the unit number
    call ftclos(unit, status)

    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)


    return
  end subroutine dump_alms_KLOAD
  !=======================================================================
  subroutine write_alms_KLOAD(filename, nalms, alms, ncl, header, nlheader, extno)
    !=======================================================================
    !     Writes alms from to binary FITS file, FKH/Apr-99
    !     ncl is the number of columns, in the output fits file,
    !     either 3 or 5 (with or without errors respectively)
    !
    !    input array (real)                   FITS file
    !     alms(:,1) = l                      )
    !     alms(:,2) = m                      )---> col 1: l*l+l+m+1
    !     alms(:,3) = real(a_lm)              ---> col 2
    !     alms(:,4) = imag(a_lm)              ---> col 3
    !     alms(:,5) = real(delta a_lm)        ---> col 4
    !     alms(:,6) = imag(delta a_lm)        ---> col 5
    !
    !=======================================================================
    
    INTEGER(I4B), INTENT(IN) :: nalms, nlheader, ncl, extno
    REAL(KMAP),        INTENT(IN), DIMENSION(0:nalms-1,1:(ncl+1)) :: alms
    CHARACTER(LEN=80), INTENT(IN), DIMENSION(1:nlheader) :: header
    CHARACTER(LEN=*),  INTENT(IN)               :: filename

    INTEGER(I4B), DIMENSION(:), allocatable :: lm
    INTEGER(I4B) ::  status,unit,blocksize,bitpix,naxis,naxes(1)
    INTEGER(I4B) ::  i,hdutype, lmax, mmax, lmin
    LOGICAL(LGT) ::  simple,extend
    CHARACTER(LEN=80) :: comment

    INTEGER(I4B), PARAMETER :: MAXDIM = MAXDIM_TOP !number of columns in the extension
    INTEGER(I4B) :: nrows, npix, tfields, varidat, repeat
    INTEGER(I4B) :: frow,  felem, colnum, stride, istart, iend, k
    CHARACTER(LEN=20) :: ttype(MAXDIM), tform(MAXDIM), tunit(MAXDIM), extname
    CHARACTER(LEN=10) ::  card
    CHARACTER(LEN=2) :: stn
    INTEGER(I4B) :: itn
    integer(I4B) :: l, m
    character(len=filenamelen) sfilename
    character(len=1) :: pform
    !-----------------------------------------------------------------------

    if (KMAP == SP) pform = 'E'
    if (KMAP == DP) pform = 'D'

    status=0
    unit = 140

    !     create the new empty FITS file
    blocksize=1

    if (extno==1) then

       call ftinit(unit,filename,blocksize,status)
       if (status > 0) call fatal_error("Error while creating file " &
            & //trim(filename) &
            & //". Check path and/or access rights.")

       !     -----------------------------------------------------
       !     initialize parameters about the FITS image
       simple=.true.
       bitpix=32     ! integer*4
       naxis=0       ! no image
       naxes(1)=0
       extend=.true. ! there is an extension

       !     ----------------------
       !     primary header
       !     ----------------------
       !     write the required header keywords
       call ftphpr(unit,simple,bitpix,naxis,naxes,0_i4b,1_i4b,extend,status)

       !     writes supplementary keywords : none

       !     write the current date
       call ftpdat(unit,status) ! format ccyy-mm-dd

       !     ----------------------
       !     image : none
       !     ----------------------

       !     ----------------------
       !     extension
       !     ----------------------

    else

       !*********************************************
       !     reopen an existing file and go to the end
       !*********************************************
       ! remove the leading '!' (if any) when reopening the same file
       sfilename = adjustl(filename)
       if (sfilename(1:1) == '!') sfilename = sfilename(2:filenamelen)
       call ftopen(unit,sfilename,1_i4b,blocksize,status)
       call ftmahd(unit,extno,hdutype,status)

    endif

    !     creates an extension
    call ftcrhd(unit, status)

    !     writes required keywords
    nrows    = nalms  ! naxis1
    tfields  = ncl
    repeat   = 1
    tform(1)(1:2)='1J' ! necessary with ifort 11.1.076
    do i=2, ncl
       tform(i)(1:2) = '1'//pform
    enddo
    ttype(1)(1:15) = 'index=l^2+l+m+1'
    ttype(2)(1:15) = 'alm (real)     '
    ttype(3)(1:15) = 'alm (imaginary)' 
    if (ncl>3) then
       ttype(4)(1:17)  = 'error (real)     '
       ttype(5)(1:17)  = 'error (imaginary)'
    endif
    tunit(1:ncl) = ''      ! optional, will not appear
    extname  = ''      ! optional, will not appear
    varidat  = 0
    call ftphbn(unit, nrows, tfields, ttype, tform, tunit, &
         &     extname, varidat, status)

    !     write the header literally, putting TFORMi at the desired place
    do i=1,nlheader
       card = header(i)
       if (card(1:5) == 'TTYPE') then ! if TTYPEi is explicitely given
          stn = card(6:8)
          read(stn,'(i3)') itn
          ! discard at their original location:
          call ftdkey(unit,'TTYPE'//stn,status)  ! old TTYPEi and  ! remove
          status = 0
          call ftdkey(unit,'TFORM'//stn,status)  !     TFORMi
          status = 0
          if (itn <= tfields) then ! only put relevant information 2008-08-27
             call putrec(unit,header(i), status)           ! write new TTYPE1
             status = 0
             if (itn==1) then
                comment = 'data format of field: 4-byte INTEGER'
             else
                if (KMAP == SP) comment = 'data format of field: 4-byte REAL'
                if (KMAP == DP) comment = 'data format of field: 8-byte REAL'
             endif
             call ftpkys(unit,'TFORM'//stn,tform(itn),comment,status) ! and write new TFORM1 right after
          endif
       elseif (header(i)/=' ') then
          call putrec(unit,header(i), status)
       endif
       status = 0
    enddo

    lmax = nint(maxval(alms(:,1)))
    lmin = nint(minval(alms(:,1)))
    mmax = nint(maxval(alms(:,2)))
    call ftukyj(unit, 'MIN-LPOL', lmin, 'Minimum L multipole order',  status)
    call ftukyj(unit, 'MAX-LPOL', lmax, 'Maximum L multipole order',  status)
    call ftukyj(unit, 'MAX-MPOL', mmax, 'Maximum M multipole degree', status)

    !     write the extension by blocks of rows ! EH, Dec 2004
    felem  = 1  ! starting position (element)
!!!    stride = 1000 ! 1000 rows at a time
    call ftgrsz(unit, stride, status) ! find optimal stride in rows
    stride = max( stride, 1)
    allocate(lm(0:stride-1))
    do k = 0, (nalms-1)/(stride * repeat)
       istart = k * (stride * repeat)
       iend   = min(nalms, istart + stride * repeat) - 1
       do i = istart, iend ! recode the (l,m) -> lm mapping, EH, 2002-08-08
          l = nint(alms(i,1))
          m = nint(alms(i,2))
          lm(i-istart) = (l+1)*l + m + 1
       enddo

       frow = istart/repeat + 1
       npix = iend - istart + 1
       call ftpclj(      unit, 1_i4b,  frow, felem, npix, lm,                         status)
       do colnum = 2, ncl
          call f90ftpcl_(unit, colnum, frow, felem, npix, alms(istart:iend,colnum+1), status)
       enddo
    enddo
    deallocate(lm)

    !     ----------------------
    !     close and exit
    !     ----------------------

    !     close the file and free the unit number
    call ftclos(unit, status)

    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)

    return
  end subroutine write_alms_KLOAD
  !=======================================================================
  subroutine read_alms_KLOAD(filename, nalms, alms, ncl, header, nlheader, extno) 
    !=======================================================================
    !     Read a FITS file
    !
    !     slightly modified to deal with vector column 
    !     in binary table       EH/IAP/Jan-98
    !
    !     Used by synfast when reading a binary file with alms for cons.real.
    !                        FKH/Apr-99
    !
    !     called by fits2alms
    !     Jan 2005, EH, improved for faster reading
    !=======================================================================
    CHARACTER(LEN=*),               INTENT(IN) :: filename
    INTEGER(I4B),                        INTENT(IN) :: nalms, ncl,nlheader,extno
    REAL(KMAP), DIMENSION(0:nalms-1,1:(ncl+1)), INTENT(OUT) :: alms
    CHARACTER(LEN=80), INTENT(OUT), DIMENSION(1:nlheader) :: header
    REAL(KMAP)                                        :: nullval
    LOGICAL(LGT)                                    ::  anynull

    INTEGER(I4B), DIMENSION(:), allocatable :: lm
    INTEGER(I4B) :: status,unit,readwrite,blocksize,naxes(2),nfound, naxis
    INTEGER(I4B) :: npix
    CHARACTER(LEN=80) :: comment ! , record
    LOGICAL(LGT) :: extend
    INTEGER(I4B) :: nmove, hdutype ! , nkeys , nspace
    INTEGER(I4B) :: frow, imap
    INTEGER(I4B) :: datacode, repeat, width
    integer(I4B) :: i, l, m
    integer(i4b) :: nrow2read, nelem
    integer(i8b) :: i0, i1

    INTEGER(I4B), PARAMETER :: MAXDIM = MAXDIM_TOP !number of columns in the extension
    INTEGER(I4B) :: nrows, tfields, varidat
    CHARACTER(LEN=20) :: ttype(MAXDIM), tform(MAXDIM), tunit(MAXDIM), extname

    !-----------------------------------------------------------------------
    status=0
    header=''
    comment=''
    ttype=''
    tform=''
    tunit=''
    extname=''
    unit = 139
    naxes(1) = 1
    naxes(2) = 1
    nfound = -1
    anynull = .false.
    alms=0.
    readwrite=0
    call ftopen(unit,filename,readwrite,blocksize,status)
    if (status > 0) then 
       call printerror(status)
       call fatal_error("Aborting.")
    endif
    !     -----------------------------------------

    !     determines the presence of image
    call ftgkyj(unit,'NAXIS', naxis, comment, status)
    if (status > 0) call printerror(status)

    !     determines the presence of an extension
    call ftgkyl(unit,'EXTEND', extend, comment, status)
    if (status > 0) status = 0 ! no extension : 
    !     to be compatible with first version of the code

    call assert (extend, 'No extension!')
    nmove = +extno
    call ftmrhd(unit, nmove, hdutype, status)
    !cc         write(*,*) hdutype

    call assert(hdutype==2, 'this is not a binary table')

    header = ""
    call get_clean_header( unit, header, filename, status)

    !        reads all the keywords
    call ftghbn(unit, MAXDIM, &
         &        nrows, tfields, ttype, tform, tunit, extname, varidat, &
         &        status)

    if (tfields<ncl) then
       print *,'found ',tfields,' columns in the file'
       print *,'expected ',ncl
       call fatal_error
    endif
    !        finds the bad data value
!     if (KMAP == SP) call ftgkye(unit,'BAD_DATA',nullval,comment,status)
!     if (KMAP == DP) call ftgkyd(unit,'BAD_DATA',nullval,comment,status)
    call f90ftgky_(unit, 'BAD_DATA', nullval, comment, status)
    if (status == 202) then ! bad_data not found
       if (KMAP == SP) nullval = s_bad_value ! default value
       if (KMAP == DP) nullval = d_bad_value ! default value
       status = 0
    endif
    !parse TFORM keyword to find out the length of the column vector
    call ftbnfm(tform(1), datacode, repeat, width, status)
    npix = nrows * repeat
!    if (npix /= nalms) then
    if (npix > nalms) then
       print *,'found ',npix,' alms'
!       print *,'expected ',nalms
       print *,'expected ',nalms,' or less'
       call fatal_error
    endif

    call ftgrsz(unit, nrow2read, status)
    nrow2read = max(nrow2read, 1)
    nelem = nrow2read * repeat
    i0 = 0_i8b
    allocate(lm(0:nelem-1))
    do frow = 1, nrows, nrow2read
       i1 = min(i0 + nrow2read * repeat, int(npix,i8b)) - 1_i8b
       nelem = i1 - i0 + 1
       ! first column -> index
       call ftgcvj(unit, 1_i4b, frow, 1_i4b, nelem, i_bad_value, &
         &        lm(0), anynull, status)
       call assert(.not. anynull, 'There are undefined values in the table!')
       ! other columns -> a(l,m)
       do imap = 2, ncl
          call f90ftgcv_(unit, imap, frow, 1_i4b, nelem, nullval, &
               &        alms(i0:i1,imap+1), anynull, status)
          call assert (.not. anynull, 'There are undefined values in the table!')
       enddo
       ! recoding of the mapping, EH, 2002-08-08
       ! lm = l*l + l + m + 1
       do i = i0, i1
          l = int(sqrt(   real(lm(i-i0)-1, kind = DP)  ) )
          m = lm(i-i0) - l*(l+1) - 1
          ! check that round-off did not screw up mapping
          if (abs(m) > l) then
             print*,'Inconsistent l^2+l+m+1 -> l,m mapping'
             print*,l, m, l*(l+1)+m+1, lm(i-i0)
             call fatal_error
          endif
          alms(i,1) = real( l, kind=KMAP)
          alms(i,2) = real( m, kind=KMAP)
       enddo
       i0 = i1 + 1_i8b
    enddo
    deallocate(lm)
    ! sanity check
    if (i0 /= npix) then
       print*,'something wrong during piece-wise reading'
       call fatal_error
    endif

    !     close the file
    call ftclos(unit, status)


    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)
    return
  end subroutine read_alms_KLOAD
  !**************************************************************************
  SUBROUTINE read_bintod_KLOAD(filename, tod, npixtot, ntods, firstpix, nullval, anynull, &
                          header, extno, units)
  !**************************************************************************
    !=======================================================================
    !     Read a FITS file
    !
    !     slightly modified to deal with vector column (ie TFORMi = '1024E')
    !     in binary table       EH/IAP/Jan-98
    !
    !     This routine is used for reading TODS by anafast.
    !     Modified to start at a given pix numb OD & RT 02/02
    !     Modified to handle huge array (npix_tot > 2^32) OD & EH 07/02
    !      2002-07-08 : bugs correction by E.H. 
    !      2007-01-31 : added units optional output
    !=======================================================================
    
    IMPLICIT NONE
    
    CHARACTER(LEN=*),               INTENT(IN)  :: filename
    INTEGER(I8B)   ,                INTENT(IN)  :: npixtot,firstpix
    INTEGER(I4B),                   INTENT(IN)  :: ntods
    REAL(KMAP), DIMENSION(0:,1:),   INTENT(OUT) :: tod
    REAL(KMAP),                     INTENT(OUT) :: nullval
    LOGICAL(LGT),                   INTENT(OUT) :: anynull
    character(len=*), dimension(1:),intent(out), optional :: header
    INTEGER(I4B),                   INTENT(IN),  OPTIONAL :: extno
    character(LEN=*), dimension(1:),intent(OUT), optional :: units

    
    integer(I4B) :: i, nl_header, len_header, nl_units, len_units
    INTEGER(I4B) :: status,unit,readwrite,blocksize,naxes(2),nfound, naxis
    INTEGER(I4B) :: npix_32 !,firstpix_32
    CHARACTER(LEN=80) :: comment
    LOGICAL(LGT) :: extend
    INTEGER(I4B) :: nmove, hdutype
    INTEGER(I4B) :: column, frow, itod
    INTEGER(I4B) :: datacode, repeat, width
    
    INTEGER(I4B), PARAMETER :: MAXDIM = MAXDIM_TOP !number of columns in the extension
    INTEGER(I4B) :: nrows, tfields, varidat,felem
    CHARACTER(LEN=20), dimension(1:MAXDIM) :: ttype, tform,tunit
    CHARACTER(LEN=20) :: extname
 
    INTEGER(I8B) :: q,iq,npix_tmp,firstpix_tmp, i0, i1
    
    !-----------------------------------------------------------------------
    status=0
    
    unit = 138
    naxes(1) = 1
    naxes(2) = 1
    nfound = -1
    anynull = .FALSE.
    
    nl_header = 0
    if (present(header)) then
       nl_header = size(header)
       len_header = 80
    endif

    nl_units = 0
    if (present(units)) then
       nl_units = size(units)
       len_units = min(80,len(units(1))) ! due to SUN compiler bug
    endif

    readwrite=0
    CALL ftopen(unit,filename,readwrite,blocksize,status)
    IF (status .GT. 0) CALL printerror(status)
    !     -----------------------------------------
    
    !     determines the presence of image
    CALL ftgkyj(unit,'NAXIS', naxis, comment, status)
    IF (status .GT. 0) CALL printerror(status)
    
    !     determines the presence of an extension
    CALL ftgkyl(unit,'EXTEND', extend, comment, status)
    IF (status .GT. 0) status = 0 ! no extension : 
    !     to be compatible with first version of the code
    
    IF (naxis .GT. 0) THEN ! there is an image
       print*,'WARNING : Image is ignored in '//trim(filename)
    ENDIF
    IF (extend) THEN ! there is an extension

       nmove = +1
       if (present(extno)) nmove = +1 + extno
       CALL ftmrhd(unit, nmove, hdutype, status)
       !cc         write(*,*) hdutype

       call assert(hdutype==2, 'this is not a binary table')

       ! reads all the keywords
       CALL ftghbn(unit, MAXDIM, &
            &        nrows, tfields, ttype, tform, tunit, extname, varidat, &
            &        status)

       IF (tfields .LT. ntods) THEN
          PRINT *,'found ',tfields,' tods in the file'
          PRINT *,'expected ',ntods
          call fatal_error
       ENDIF

       if (nl_header > 0) then
          header = ""
          status = 0
          call get_clean_header(unit, header, filename, status)
       endif

       if (nl_units > 0) then
          do i=1,nl_units
             units(i)(1:len_units) = 'unknown' ! default
          enddo
          do itod = 1, min(ntods, nl_units)
             units(itod)(1:len_units) = adjustl(tunit(itod))
          enddo
       endif

       ! finds the bad data value
!        if (KMAP == SP) CALL ftgkye(unit,'BAD_DATA',nullval,comment,status)
!        if (KMAP == DP) CALL ftgkyd(unit,'BAD_DATA',nullval,comment,status)
       call f90ftgky_(unit, 'BAD_DATA', nullval, comment, status)
       IF (status .EQ. 202) THEN ! bad_data not found
          if (KMAP == SP) nullval = s_bad_value ! default value
          if (KMAP == DP) nullval = d_bad_value ! default value
          status = 0
       ENDIF

       IF (npixtot .LT. nchunk_max) THEN

          DO itod = 1, ntods

             !parse TFORM keyword to find out the length of the column vector (repeat)
             CALL ftbnfm(tform(itod), datacode, repeat, width, status)
             frow = (firstpix)/repeat+1          ! 1 based 
             felem = firstpix-(frow-1)*repeat+1  ! 1 based 
             npix_32 = npixtot 

             !reads the columns
             column = itod
             CALL f90ftgcv_(unit, column, frow, felem, npix_32, nullval, &
                  &        tod(0:npix_32-1,itod), anynull, status)
          END DO

       ELSE

          q = (npixtot-1)/nchunk_max
          DO iq = 0,q
             IF (iq .LT. q) THEN
                npix_tmp = nchunk_max
             ELSE
                npix_tmp = npixtot - iq*nchunk_max
             ENDIF
             firstpix_tmp = firstpix + iq*nchunk_max
             npix_32 = npix_tmp
             i0 = firstpix_tmp-firstpix
             i1 = i0 + npix_tmp - 1_i8b

             DO itod = 1, ntods
                ! parse TFORM keyword to find out the length of the column vector
                CALL ftbnfm(tform(itod), datacode, repeat, width, status)
                frow = (firstpix_tmp)/repeat+1          ! 1 based 
                felem = firstpix_tmp-(frow-1)*repeat+1  ! 1 based 
                CALL f90ftgcv_(unit, itod, frow, felem, npix_32, nullval, &
                     &      tod(i0:i1,itod), anynull, status)
             END DO

          ENDDO

       ENDIF

    ELSE ! no image no extension, you are dead, man
       call fatal_error(' No image, no extension')
    ENDIF

    ! close the file
    CALL ftclos(unit, status)

    ! check for any error, and if so print out error messages
    IF (status .GT. 0) CALL printerror(status)

    RETURN

  END SUBROUTINE read_bintod_KLOAD
  !=======================================================================
  
  !======================================================================================
  SUBROUTINE write_bintabh_KLOAD(tod, npix, ntod, header, nlheader, filename, extno, firstpix, repeat)
    !======================================================================================

    ! =================================================================================
    !     Create a FITS file containing a binary table extension in the first extension
    !
    !     Designed to deal with Huge file, (n_elements > 2^31)
    !
    !     OPTIONAL NEW PARAMETERS:
    !     firstpix : position in the file of the first element to be written (starts at 0) 
    !                default value =0
    !                8-byte integer
    !                if NE 0 then suppose that the file already exists
    !
    !     repeat   : length of vector per unit rows and columns of the first binary extension
    !                default value = 12000 (\equiv 1 mn of PLANCK/HFI data)
    !                4-byte integer
    ! 
    !     OTHER PARAMETERS
    !     unchanged with respect to the standard write_bintab of the HEALPIX package except 
    !     npix which is an 8-byte integer
    !
    !     Adapted from write_bintab
    !                                           E.H. & O.D. @ IAP 07/02
    !
    !     Requires a compilation in 64 bits of the CFITSIO 
    !     Note that the flag -D_FILE_OFFSETS_BITS=64 has to be added 
    !         (cf page CFITIO 2.2 User's guide  Chap 4, section 4-13)
    ! 
    ! 2002-07-08 : bugs correction by E.H. 
    !    (uniform use of firstpix_tmp, introduction of firstpix_chunk)
    ! 2015-04-28: better handling of repeat (suggested by R. Keskitalo)
    !==========================================================================================

    USE healpix_types
    IMPLICIT NONE

    INTEGER(I8B)     , INTENT(IN)           :: npix
    INTEGER(I8B)     , INTENT(IN), OPTIONAL :: firstpix
    INTEGER(I4B)     , INTENT(IN), OPTIONAL :: repeat
    INTEGER(I4B)     , INTENT(IN)           :: ntod,nlheader
    REAL(KMAP)       , INTENT(IN), DIMENSION(0:npix-1,1:ntod) :: tod
    CHARACTER(LEN=80), INTENT(IN), DIMENSION(1:nlheader)      :: header
    CHARACTER(LEN=*),  INTENT(IN)           :: filename
    INTEGER(I4B)    , INTENT(IN)     , OPTIONAL :: extno

    INTEGER(I4B) :: status,unit,blocksize,bitpix,naxis,naxes(1),repeat_fits
    INTEGER(I4B) :: i,npix_32
    LOGICAL(LGT) :: simple,extend
    CHARACTER(LEN=80) :: comment, ch
    integer(I8B) :: repeat_tmp

    INTEGER(I4B), PARAMETER :: MAXDIM = MAXDIM_TOP !number of columns in the extension
    INTEGER(I4B)      :: nrows,tfields,varidat
    INTEGER(I4B)      :: frow,felem,colnum,readwrite,width,datacode,hdutype
    CHARACTER(LEN=20) :: ttype(MAXDIM), tform(MAXDIM), tunit(MAXDIM), extname
    CHARACTER(LEN=10) :: card
    CHARACTER(LEN=2)  :: stn
    INTEGER(I4B)      :: itn  

    real(KMAP), dimension(:), allocatable :: padding
    integer(I8B) :: lastpix
    integer(I4B) :: lengap

    INTEGER(I4B)      :: extno_i
    character(len=filenamelen) :: sfilename
    INTEGER(I8B) :: q,iq,npix_tmp,firstpix_tmp,firstpix_chunk, i0, i1
    character(len=1) :: pform
    character(len=*), parameter :: code='write_bintabh'
    !-----------------------------------------------------------------------

    if (KMAP == SP) pform = 'E'
    if (KMAP == DP) pform = 'D'

    IF (.NOT. PRESENT(repeat) ) THEN 
       repeat_tmp = 1
       if (mod(npix,1024_i8b) == 0) then
          repeat_tmp = 1024
       elseif (npix >= 12000) then
          repeat_tmp = 12000 
       endif
    ELSE 
       repeat_tmp = repeat
    ENDIF
    IF (.NOT. PRESENT(firstpix) ) THEN 
       firstpix_tmp = 0 
    ELSE 
       firstpix_tmp = firstpix
    ENDIF

    extno_i = 0
    if (present(extno)) extno_i = extno

    status=0
    unit = 137
    blocksize=1

    ! remove the leading '!' (if any) when reopening the same file
    sfilename = adjustl(filename)
    if (sfilename(1:1) == '!') sfilename = sfilename(2:filenamelen)

    ! create the new empty FITS file

    IF (firstpix_tmp .EQ. 0) THEN

       if (extno_i == 0) then 
          CALL ftinit(unit,filename,blocksize,status)
          if (status > 0) call fatal_error("Error while creating file " &
               & //trim(filename) &
               & //". Check path and/or access rights.")

          ! -----------------------------------------------------
          ! Initialize parameters about the FITS image
          simple=.TRUE.
          bitpix=32     ! integer*4
          naxis=0       ! no image
          naxes(1)=0
          extend=.TRUE. ! there is an extension

          !     ----------------------
          !     primary header
          !     ----------------------
          !     write the required header keywords
          CALL ftphpr(unit,simple,bitpix,naxis,naxes,0_i4b,1_i4b,extend,status)

          !     writes supplementary keywords : none

          !     write the current date
          CALL ftpdat(unit,status) ! format ccyy-mm-dd

       !     ----------------------
       !     image : none
       !     ----------------------

       !     ----------------------
       !     extension
       !     ----------------------
       else

          !*********************************************
          !     reopen an existing file and go to the end
          !*********************************************
          call ftopen(unit,sfilename,1_i4b,blocksize,status)
          call ftmahd(unit,1_i4b+extno_i,hdutype,status)

       endif

       !     creates an extension
       CALL ftcrhd(unit, status)

       !     writes required keywords
       ! if (npix < repeat_tmp) repeat_tmp = npix ! 2015-04-28

       ! nrows    = npix / repeat_tmp + 1 ! naxis1
       nrows    = (npix - 1_i8b) / repeat_tmp + 1_i4b ! 2015-04-28
       tfields  = ntod
       WRITE(ch,'(i8)') repeat_tmp
!       tform(1:ntod)  = TRIM(ADJUSTL(ch))//pform ! does not work with Ifort, EH, 2006-04-04
       ch = TRIM(ADJUSTL(ch))//pform
       tform(1:ntod)  = ch

!        IF (npix .LT. repeat_tmp) THEN  ! 2015-04-28
!           nrows = npix
!           ch = '1'//pform
!           tform(1:ntod) = ch
!        ENDIF
       ttype(1:ntod) = 'simulation'   ! will be updated
       tunit(1:ntod) = ''      ! optional, will not appear

       extname  = ''      ! optional, will not appear
       varidat  = 0

       CALL ftphbn(unit, nrows, tfields, ttype, tform, tunit, &
            &     extname, varidat, status)

       !     write the header literally, putting TFORM1 at the desired place
       if (KMAP == SP) comment = 'data format of field: 4-byte REAL'
       if (KMAP == DP) comment = 'data format of field: 8-byte REAL'
       DO i=1,nlheader
          card = header(i)
          IF (card(1:5) == 'TTYPE') THEN ! if TTYPE1 is explicitely given
             stn = card(6:8)
             READ(stn,'(i3)') itn
             ! discard at their original location:
             call ftdkey(unit,'TTYPE'//stn,status)  ! old TTYPEi and
             status = 0
             call ftdkey(unit,'TFORM'//stn,status)  !     TFORMi
             status = 0
             if (itn <= tfields) then ! only put relevant information 2008-08-27
                call putrec(unit,header(i), status)    ! write new TTYPE1
                status = 0
                call ftpkys(unit,'TFORM'//stn,tform(itn),comment,status) ! and write new TFORM1 right after
                !CALL ftmcrd(unit,'TTYPE'//stn,'COMMENT',status)  ! old TTYPEi and 
                !CALL ftmcrd(unit,'TFORM'//stn,'COMMENT',status)  !     TFORMi
                !CALL ftprec(unit,header(i), status)           ! write new TTYPE1
                !CALL ftpkys(unit,'TFORM'//stn,tform(1),comment,status) ! and write new TFORM1 right after
             endif
          ELSEIF (header(i).NE.' ') THEN
             call putrec(unit,header(i), status)
!              CALL ftprec(unit,header(i), status)
          ENDIF
          status = 0
       ENDDO

    ELSE
       ! The file already exists
       readwrite=1
       CALL ftopen(unit,sfilename,readwrite,blocksize,status)
       CALL ftmahd(unit,2_i4b+extno_i,hdutype,status) 

       CALL ftgkys(unit,'TFORM1',tform(1),comment,status)
       CALL ftbnfm(tform(1),datacode,repeat_fits,width,status)

       IF (repeat_tmp .NE. repeat_fits) THEN
          if (present(repeat)) then 
             write(*,'(a)') code//'> WARNING: In file '//trim(sfilename)
             write(*,'(a)') &
                  & code//'>  user provided REPEAT value (' // &
                  & trim(adjustl(string(repeat_tmp)))       // &
                  & ') differs from value read from file (' // &
                  & trim(adjustl(string(repeat_fits)))      // &  
                  & ').'
             write(*,'(a)') code//'> The latter will be used.'
!           else
!              write(*,'(a,i0.0,a)') &
!                   & code//'> WARNING: REPEAT value read from file (', &
!                   & repeat_fits, &
!                   & ') will be used.'
          endif
          repeat_tmp = repeat_fits
       ENDIF

    ENDIF


    IF (npix .LT. nchunk_max) THEN ! data is small enough to be written in one chunk

       frow = (firstpix_tmp)/repeat_tmp + 1
       felem = firstpix_tmp-(frow-1)*repeat_tmp + 1
       npix_32 = npix 

       DO colnum = 1, ntod
          call f90ftpcl_(unit, colnum, frow, felem, npix_32, tod(0:npix_32-1,colnum), status)
       END DO

    ELSE ! data has to be written in several chunks

       q = (npix-1)/nchunk_max
       DO iq = 0,q
          IF (iq .LT. q) THEN
             npix_tmp = nchunk_max
          ELSE
             npix_tmp = npix - iq*nchunk_max
          ENDIF
          i0 = iq*nchunk_max
          i1 = i0 + npix_tmp - 1_i8b
          firstpix_chunk = firstpix_tmp + i0
          frow  = (firstpix_chunk)/repeat_tmp + 1
          felem =  firstpix_chunk-(frow-1)*repeat_tmp + 1
          npix_32 = npix_tmp
          DO colnum = 1, ntod
             call f90ftpcl_(unit, colnum, frow, felem, npix_32, &
                  &          tod(i0:i1, colnum), status)
          END DO
       ENDDO

    ENDIF

    ! pad entry if necessary ! 2015-04-28
    lastpix = firstpix_tmp + npix  ! number of pixels written above
    lengap = modulo(lastpix, repeat_tmp) ! remaining gap
    if (lengap > 0) then
       firstpix_tmp = lastpix
       npix_32 = repeat_tmp - lengap
       frow    = (firstpix_tmp)/repeat_tmp + 1
       felem   =  firstpix_tmp-(frow-1)*repeat_tmp + 1
       allocate(padding(0:npix_32-1))
       if (KMAP == SP) padding(:) = HPX_SBADVAL
       if (KMAP == DP) padding(:) = HPX_DBADVAL
       do colnum = 1, ntod
          call f90ftpcl_(unit, colnum, frow, felem, npix_32, padding, status)
       enddo
       deallocate(padding)
    endif

    ! ----------------------
    ! close and exit
    ! ----------------------

    ! close the file and free the unit number
    CALL ftclos(unit, status)

    ! check for any error, and if so print out error messages
    IF (status .GT. 0) CALL printerror(status)

    RETURN

  END SUBROUTINE write_bintabh_KLOAD
  ! ==============================================================================

!******************************************************************************
!---------------------------------------
!  reads a FITS file containing a list of 
!  ring-based or pixel-based quadrature weights
!  and turns it into a (ring-ordered) full sky map
!
! adapted from unfold_weights.pro
! 2018-05-25
!---------------------------------------
  subroutine unfold_weightsfile_KLOAD(w8file, w8map)

    USE pix_tools,      only: nside2npix, nside2npweights

    character(len=*), intent(in)                  :: w8file
    real(KMAP),       intent(out), dimension(0:)  :: w8map
    real(KMAP),       allocatable, dimension(:,:) :: w8list
    integer(i4b)                :: nside, status, won
    integer(i8b)                :: nfw8, npw8, nrw8
    character(len=*), parameter :: code = 'unfold_weightsfile'

    nfw8 = getsize_fits(w8file, nside=nside)
    npw8 = nside2npweights(nside)
    nrw8 = 2*nside

    won = 0
    if (nfw8 == nrw8) won = 1
    if (nfw8 == npw8) won = 2
    if (won == 0) then
       print*,'Weights file = '//trim(w8file)
       print*,'contains ',nfw8,' weights'
       print*,'for Nside    = ',nside
       print*,'Expected either ',nrw8,' or ',npw8
       call fatal_error
    endif

    allocate(w8list(0:nfw8-1,1:1), stat=status)
    call assert_alloc(status,code,"w8list")
    call input_map(w8file, w8list, nfw8, 1)
    call unfold_weightslist(nside, won, w8list(:,1), w8map)
    deallocate(w8list)
    
    return
  end subroutine unfold_weightsfile_KLOAD
!---------------------------------------
!  turns a list of ring-based or pixel-based quadrature weights 
!  into a (ring-ordered) full sky map
!
! adapted from unfold_weights.pro
! 2018-05-18
! 2018-09-12: promoted qpix to I8B to avoid consistency problems during compilation with g95  
!---------------------------------------
  subroutine unfold_weightslist_KLOAD(nside, won, w8list, w8map)

    USE pix_tools,      only: nside2npix, nside2npweights
    USE long_intrinsic, only: long_size

    integer(i4b),            intent(in)  :: nside, won
    real(KMAP), dimension(0:), intent(in)  :: w8list
    real(KMAP), dimension(0:), intent(out) :: w8map

    integer(i8b) :: np, npix
    integer(i8b) :: nf, nw8, pnorth, psouth, vpix, qp4, p, qpix
    integer(i4b) :: ring, odd, shifted, rpix, wpix, j4, n4, it

    ! test nside
    npix = nside2npix(nside)
    if (npix <= 0) then 
       print*,'Unvalid Nside = ',nside
       call fatal_error
    endif

    ! test won
    if (won < 1 .or. won > 2) then
       print*,'Expected either 1 or 2, got ',won
       call fatal_error
    endif

    ! test input list size
    nf  = long_size(w8list)
    if (won == 1) then
       nw8 = 2 * nside
    else
       nw8 = nside2npweights(nside)
    endif
    if (nf /= nw8) then
       print*,'Expected in input weight list',nw8
       print*,'got ',nf
       call fatal_error
    endif

    ! test output map size
    np = long_size(w8map)
    if (np /= npix) then
       print*,'Expected in output weight map',npix
       print*,'got ',np
       call fatal_error
    endif

    ! do actual unfolding
    if (won == 1) then

       pnorth = 0_i8b
       n4 = 4 * nside
       do ring=1,n4-1                ! loop on all rings
          it = min(ring, n4 - ring)  ! ring index, counting from closest pole
          qp4= 4_i8b * min(it, nside)    ! pixels / ring
          w8map(pnorth:pnorth+qp4-1) = w8list(it-1)
          pnorth = pnorth + qp4
       enddo

    else
       
       pnorth = 0_i8b            ! position in expanded weights
       vpix   = 0_i8b            ! position in compress list
       do ring=0, 2*nside-1
           qpix    = min(ring+1, nside) ! 1/4 of number of pixels per ring
           shifted = 0
           if (ring < (nside-1) .or. iand(ring+nside,1_i4b) == 1_i4b) shifted = 1
           odd     = iand(qpix,1_i8b)
           qp4     = 4_i8b*qpix       ! number of pixels per ring
           ! fill the weight map
           do p=0,qp4-1
              j4 = mod(p, qpix) ! g95: p and qpix must be of same type
              rpix = min(j4, qpix - shifted - j4) ! seesaw
              w8map(pnorth+p) = w8list(vpix + rpix)
           enddo
       
           if (ring < 2*nside-1) then  ! south part (except equator)
               psouth     = npix - pnorth - qp4
               w8map(psouth:psouth+qp4-1) = w8map(pnorth:pnorth+qp4-1)
           endif
           ! locations on next ring
           pnorth  = pnorth + qp4
           wpix    = (qpix+1)/2  + 1 - ior(odd, shifted)
           vpix    = vpix + wpix
       enddo
    endif

    ! add 1 to get final weight
    w8map = w8map + 1.0_KMAP

    return
  end subroutine unfold_weightslist_KLOAD

  !=======================================================================
  ! read_fits_partial(filename, pixel, cutmap, [header, extno])
  !   routine to read FITS file with cut sky data : PIXEL, ?, ?, ?, ?, ...
  !=======================================================================
#ifndef NO64BITS
  !=======================================================================
  subroutine read_fits_partial4_KLOAD(filename, pixel, cutmap, header, extno)
    !=======================================================================
    use long_intrinsic, only: long_size
    character(len=*),                     intent(in)            :: filename
    integer(I4B),     dimension(0:),      intent(out)           :: pixel
    real(KMAP),       dimension(0:,1:),   intent(out)           :: cutmap
    character(len=*), dimension(1:),      intent(out), optional :: header
    integer(I4B),                         intent(in),  optional :: extno
    !
    integer(i8b),     dimension(:), allocatable :: pixel8
    integer(i8b) :: npix
    character(len=*), parameter :: code = 'read_fits_partial'
    !---------------------------------------------------------------------
    npix = long_size(pixel)
    if (npix > MAX_I4b) then
       call fatal_error(code//': PIXEL must be of type I8B to read '//trim(filename))
    endif
    allocate(pixel8(0:npix-1))
    call read_fits_partial8_KLOAD(filename, pixel8, cutmap, header, extno)
    if (maxval(pixel8) > MAX_I4B) then
       call fatal_error(code//': PIXEL must be of type I8B to read '//trim(filename))       
    endif
    pixel(0:npix-1) = pixel8(0:npix-1)
    deallocate(pixel8)
    return

  end subroutine read_fits_partial4_KLOAD
#endif
  !=======================================================================
  subroutine read_fits_partial8_KLOAD(filename, pixel, cutmap, header, extno)
    !=======================================================================
    use long_intrinsic, only: long_size
    character(len=*),                     intent(in)            :: filename
    integer(I8B),     dimension(0:),      intent(out)           :: pixel
    real(KMAP),       dimension(0:,1:),   intent(out)           :: cutmap
    character(len=*), dimension(1:),      intent(out), optional :: header
    integer(I4B),                         intent(in),  optional :: extno

    integer(I8B) :: obs_npix, npixtot, npix

    integer(I4B), parameter :: MAXDIM = MAXDIM_TOP !number of columns in the extension
    integer(I4B) :: blocksize, datacode
    integer(I4B) :: firstpix, frow, hdutype
    integer(I4B) :: naxis, nfound, nmove, nrows, nmaps, nmr, i, npix_4
    integer(I4B) :: readwrite, status, tfields, unit, varidat, width
    integer(I4B) :: repeat, repeat1, repeat2
    logical(LGT) :: anynull, extend
    character(len=20), dimension(MAXDIM) :: ttype, tform, tunit
    character(len=20) :: extname
    character(len=80) :: comment
    real(KMAP) :: nullval
    !=======================================================================

    npixtot = min(long_size(pixel),long_size(cutmap,1))
    nmaps   = size(cutmap,2)
    status=0

    unit = 150
    nfound = -1
    anynull = .false.

    readwrite=0
    call ftopen(unit,filename,readwrite,blocksize,status)
    if (status > 0) call printerror(status)

    !     determines the presence of image
    call ftgkyj(unit,'NAXIS', naxis, comment, status)
    if (status > 0) call printerror(status)
    if (naxis > 0) then ! there is an image
       print*,'an image was found in the FITS file '//filename
       print*,'... it is ignored.'
    endif

    !     determines the presence of an extension
    call ftgkyl(unit,'EXTEND', extend, comment, status)
    if (status > 0) then 
       print*,'extension expected and not found in FITS file '//filename
       print*,'abort code'
       call fatal_error
    endif

    nmove = +1
    if (present(extno)) nmove = +1 + extno
    call ftmrhd(unit, nmove, hdutype, status)

    call assert (hdutype==2, 'this is not a binary table')

    ! reads all the FITS related keywords
    call ftghbn(unit, MAXDIM, &
         &        nrows, tfields, ttype, tform, tunit, extname, varidat, &
         &        status)

    if (.not.  (trim(ttype(1)) == 'PIXEL' &
         & .or. trim(ttype(1)) == 'PIX'  ) ) call fatal_error('did not find PIXEL column in '//filename)
    

    if (present(header)) then 
       header = ""
       status = 0
       !call fitstools_mp_get_clean_header(unit, header, filename, status)
       call get_clean_header(unit, header, filename, status)
    endif

!     if (present(units)) then 
!        ! the second column contains the SIGNAL, for which we 
!        ! already read the units from the header
!        units = adjustl(tunit(2))
!     endif

    !        finds the bad data value
    call f90ftgky_(unit, 'BAD_DATA', nullval, comment, status)
    if (status == 202) then ! bad_data not found
       if (KMAP == SP) nullval = s_bad_value ! default value
       if (KMAP == DP) nullval = d_bad_value ! default value
       status = 0
    endif

    frow = 1
    firstpix = 1
    !        parse TFORM keyword to find out the length of the column vector
    repeat1 = 1
    repeat2 = 1
    call ftbnfm(tform(1), datacode, repeat1, width, status)
    if (tfields > 1) call ftbnfm(tform(2), datacode, repeat2, width, status)
    repeat = max(repeat1, repeat2)

    !call ftgkyj(unit,'OBS_NPIX',obs_npix,comment,status) ! I4B
    !call ftgkyk(unit,'OBS_NPIX',obs_npix,comment,status) ! I8B
    call f90ftgky_(unit,'OBS_NPIX',obs_npix,comment,status) ! generic
    if (status == 202) then ! obs_npix not found
       obs_npix = int(nrows, kind=i8b) * repeat
       status = 0
       print*,'Keyword OBS_NPIX not found in '//trim(filename)
       print*,'Will assume OBS_NPIX='//trim(string(obs_npix,format='(i0.0)'))//' instead.' ! necessary to do a printout with gfortran (on M1), 2022-07-06
    endif

    npix   = min(npixtot, obs_npix)
    npix_4 = npix
    nmr    = min(nmaps,   tfields-1)
    !call ftgcvj   (unit, 1_i4b, frow, firstpix, npix_4, i_bad_value, pixel(0), anynull, status) ! I4B
    !call ftgcvk   (unit, 1_i4b, frow, firstpix, npix_4, l_bad_value, pixel(0), anynull, status) ! I8B
    call f90ftgcv_(unit, 1_i4b, frow, firstpix, npix_4, l_bad_value, pixel, anynull, status) ! generic
    do i=1,nmr
       call f90ftgcv_(unit, i+1_i4b, frow, firstpix, npix_4, nullval, cutmap(0:npix-1,i), anynull, status)
    enddo

    !     close the file
    call ftclos(unit, status)

    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)

    return
  end subroutine read_fits_partial8_KLOAD


  !=======================================================================
  ! writes a fits file for (polarized) sky data set with columns:
  !  PIXEL, TEMPERATURE, [Q_POLARISATION, U_POLARISATION]
  ! 
  ! write_fits_partial(filename, pixel, cutmap, &
  !         [, header, coord, nside, order, units, extno])
  !=======================================================================
  subroutine write_fits_partial4_KLOAD(filename, pixel, cutmap, &
       &                     header, coord, nside, order, units, extno)
    !=======================================================================
    character(len=*),                   intent(in)           :: filename
    integer(I4B),     dimension(0:),    intent(in)           :: pixel
    real(KMAP),       dimension(0:,1:), intent(in)           :: cutmap
    character(len=*), dimension(1:),    intent(in), optional :: header
    integer(I4B),                       intent(in), optional :: nside, order
    character(len=*),                   intent(in), optional :: coord, units
    integer(I4B),                       intent(in), optional :: extno !, polarisation
    character(len=*), parameter :: routine = 'write_fits_partial'
    call sub_write_fits_partial_KLOAD(filename, cutmap, pixel4=pixel, &
         header=header, coord=coord, nside=nside, order=order, units=units, extno=extno)
  end subroutine write_fits_partial4_KLOAD
  !=======================================================================
#ifndef NO64BITS
  subroutine write_fits_partial8_KLOAD(filename, pixel, cutmap, &
       &                     header, coord, nside, order, units, extno)
    !=======================================================================
    character(len=*),                   intent(in)           :: filename
    integer(I8B),     dimension(0:),    intent(in)           :: pixel
    real(KMAP),       dimension(0:,1:), intent(in)           :: cutmap
    character(len=*), dimension(1:),    intent(in), optional :: header
    integer(I4B),                       intent(in), optional :: nside, order
    character(len=*),                   intent(in), optional :: coord, units
    integer(I4B),                       intent(in), optional :: extno !, polarisation
    character(len=*), parameter :: routine = 'write_fits_partial'
    call sub_write_fits_partial_KLOAD(filename, cutmap, pixel8=pixel, &
         header=header, coord=coord, nside=nside, order=order, units=units, extno=extno)
  end subroutine write_fits_partial8_KLOAD
#endif
  !=======================================================================
  subroutine sub_write_fits_partial_KLOAD(filename, cutmap, &
       &                     pixel4, pixel8, &
       &                     header, coord, nside, order, units, extno)
    !=======================================================================
    use long_intrinsic, only: long_size
    use pix_tools,      only: nside2npix
    use misc_utils,     only: assert

    character(len=*),                   intent(in)           :: filename
    real(KMAP),       dimension(0:,1:), intent(in)           :: cutmap
    integer(I4B),     dimension(0:),    intent(in), optional :: pixel4
    integer(I8B),     dimension(0:),    intent(in), optional :: pixel8
    character(len=*), dimension(1:),    intent(in), optional :: header
    integer(I4B),                       intent(in), optional :: nside, order
    character(len=*),                   intent(in), optional :: coord, units
    integer(I4B),                       intent(in), optional :: extno !, polarisation

    character(len=*), parameter :: routine = 'write_fits_partial'
    ! --- healpix/fits related variables ----
    integer(I4B)     :: ncol, grain, nd2
    integer(I4B)     :: npix_hd, nside_hd, i
    integer(I4B)     :: maxpix, minpix
    character(len=1) :: char1, coord_usr
!    character(len=8) :: char8
    character(len=20) :: units_usr
    logical(LGT)     :: done_nside, done_order, done_coord, done_polar, polar_flag
    integer(I4B)     :: nlheader, extno_i

    ! --- cfitsio related variables ----
    integer(I4B) ::  status, unit, blocksize, bitpix, naxis, naxes(1)
    logical(LGT) ::  simple, extend
    character(LEN=80) :: svalue, comment

    integer(I4B), parameter :: MAXDIM = MAXDIM_TOP !number of columns in the extension
    integer(I4B) :: nrows, tfields, varidat
    integer(I4B) :: frow,  felem, repeat, repeatg, hdutype
    character(len=20) :: ttype(MAXDIM), tform(MAXDIM), tunit(MAXDIM), extname
    character(len=20), dimension(1:3) :: extnames
    character(len=80) ::  card
    character(len=4) :: srepeat, srepeatg
    character(len=1) :: pform1, pform2
    character(len=filenamelen) sfilename
    integer(I8B)     :: obs_npix8, npix_final, obs_npp
    integer(I4B)     :: obs_npix4, kpix, nside_final

    integer(I4B), save :: nside_old
    !=======================================================================
    if (KMAP == SP)  pform2='E'
    if (KMAP == DP)  pform2='D'
    if (present(pixel4) .eqv. present(pixel8)) then
       call fatal_error(routine//': choose either PIXEL4 or PIXEL8')
    endif
    if (present(pixel4)) then
       kpix = I4B
       pform1='J'
       obs_npp = long_size(pixel4)
    elseif (present(pixel8)) then
       kpix = I8B
       pform1='K'
       obs_npp = long_size(pixel8)
    endif
    obs_npix8 = long_size(cutmap, 1)
    nd2       =      size(cutmap, 2)
    obs_npix4 = int(obs_npix8, kind=I4B)
    call assert( obs_npix8 == obs_npp, routine//': mismatched size for PIXEL and DATA')
    ncol = 1 + nd2
    grain = 1
    status=0
    unit = 149
    blocksize=1

    nlheader = 0
    if (present(header)) nlheader = size(header)
    units_usr = ' '
    if (present(units)) units_usr = units
    extno_i = 0
    if (present(extno)) extno_i = extno
    polar_flag = .false.
    done_polar = .false.
!     if (present(polarisation)) then
!        polar_flag = (polarisation /= 0)
!        done_polar = .true. ! will ignore the header provided value of POLAR
!     endif

    if (extno_i == 0) then
       !*************************************
       !     create the new empty FITS file
       !*************************************
       call ftinit(unit,filename,blocksize,status)
       if (status > 0) call fatal_error("Error while creating file " &
            & //trim(filename) &
            & //". Check path and/or access rights.")

       !     -----------------------------------------------------
       !     initialize parameters about the FITS image
       simple=.true.
       bitpix=32     ! integer*4
       naxis=0       ! no image
       naxes(1)=0
       extend=.true. ! there is an extension

       !     ----------------------
       !     primary header
       !     ----------------------
       !     write the required header keywords
       call ftphpr(unit,simple,bitpix,naxis,naxes,0_i4b,1_i4b,extend,status)

       !     write the current date
       call ftpdat(unit,status) ! format ccyy-mm-dd

    else
       !*********************************************
       !     reopen an existing file and go to the end
       !*********************************************
       ! remove the leading '!' (if any) when reopening the same file
       sfilename = adjustl(filename)
       if (sfilename(1:1) == '!') sfilename = sfilename(2:filenamelen)
       call ftopen(unit,sfilename,1_i4b,blocksize,status)
       ! move to last extension written
       call ftmahd(unit, 1_i4b+ extno_i, hdutype, status)

    endif

       
    !     ----------------------
    !     new extension
    !     ----------------------
    call ftcrhd(unit, status)

       !     writes required keywords
       !     repeat = 1024
    repeat = 1
    if (obs_npix8 < repeat) repeat = 1
    nrows    = (obs_npix8 + repeat - 1)/ repeat ! naxis1
    write(srepeat,'(i4)') repeat
    srepeat = adjustl(srepeat)

    repeatg = repeat * grain
    write(srepeatg,'(i4)') repeatg
    srepeatg = adjustl(srepeatg)

    tfields  = ncol
    ttype(1) = 'PIXEL'
    if (nd2  == 1) then
       done_polar = .true.
       polar_flag = .false.
       ttype(2)   = 'TEMPERATURE  '
    elseif (nd2  == 3) then
       done_polar = .true.
       polar_flag = .true.
       ttype(2:4) = (/'TEMPERATURE   ' , & 
            &         'Q_POLARISATION' , &
            &         'U_POLARISATION' /)
    else
       do i=2,ncol
          ttype(i) = 'C'//string(i-1, format='(i2.2)') ! C01, C02, ...
       enddo
    endif

    tform(1)      = trim(srepeat) //pform1
    tform(2:ncol) = trim(srepeatg)//pform2

    tunit =  ' '      ! optional, will not appear
    tunit(2:ncol) = units_usr

    extname  = 'SKY_OBSERVATION'      ! default, will be overridden by user provided one if any
    ! if (polar_flag) extname = extnames(1+extno_i)
    varidat  = 0
    call ftphbn(unit, nrows, tfields, ttype, tform, tunit, &
         &     extname, varidat, status)

    call ftpcom(unit,'------------------------------------------',status)
    call ftpcom(unit,'          Pixelisation Specific Keywords    ',status)
    call ftpcom(unit,'------------------------------------------',status)
    call ftpkys(unit,'PIXTYPE','HEALPIX ',' HEALPIX Pixelisation',status)      
    call ftpkyu(unit,'NSIDE',   ' ',status) ! place holder, will be updated later on
    call ftpkyu(unit,'ORDERING',' ',status)
    call ftpkys(unit,'COORDSYS',' ',' ',status)
    call ftpcom(unit,'  G = Galactic, E = ecliptic, C = celestial = equatorial', status)   
    call ftpkyl(unit,'POLAR',polar_flag,'Polarisation included in file (T/F)',status)
    call ftpcom(unit,'------------------------------------------',status)
    call ftpcom(unit,'          Data Specific Keywords    ',status)
    call ftpcom(unit,'------------------------------------------',status)
    call ftpkys(unit,'INDXSCHM','EXPLICIT',' Indexing : IMPLICIT or EXPLICIT', status)
!     call ftpkyj(unit,'GRAIN',  grain,     ' Grain of pixel indexing',status)
!     call ftpcom(unit,'GRAIN=0 : no indexing of pixel data                         (IMPLICIT)',status)
!     call ftpcom(unit,'GRAIN=1 : 1 pixel index -> 1 pixel data                     (EXPLICIT)',status)
!     call ftpcom(unit,'GRAIN>1 : 1 pixel index -> data of GRAIN consecutive pixels (EXPLICIT)',status)
    call ftpkys(unit,'OBJECT','PARTIAL ',' Sky coverage represented by data',status)   
    if (KPIX == I4B) then
       call ftpkyj(unit,'OBS_NPIX',obs_npix4, ' Number of pixels observed and recorded',status)
    else
       call ftpkyk(unit,'OBS_NPIX',obs_npix8, ' Number of pixels observed and recorded',status)
    endif

    ! add required Healpix keywords (NSIDE, ORDER) if provided by user
    done_order = .false.
    if (present(order)) then
       !        char8 = order
       !        call ftupch(char8)
       !        if (char8(1:4) == 'RING') then
       if (order == 1) then
          call ftukys(unit, 'ORDERING','RING',' Pixel ordering scheme, either RING or NESTED',status)
          done_order = .true.
          !        elseif (char8(1:4) == 'NEST') then
       elseif (order == 2) then
          call ftukys(unit, 'ORDERING','NESTED',' Pixel ordering scheme, either RING or NESTED ',status)
          done_order = .true.
       else
          print*,'Invalid ORDER given : ',order, ' instead of 1 (RING) or 2 (NESTED)'
       endif
    endif

    done_nside = .false.
    if (present(nside)) then
       if (nside2npix(nside) > 0) then ! valid nside
          call ftukyj(unit,'NSIDE',nside,' Resolution parameter for HEALPIX', status)
          done_nside = .true.
          nside_final = nside
       else
          print*,'Invalid NSIDE given : ',nside
       endif
    endif

    ! add non required Healpix keyword (COORD)
    done_coord = .false.
    if (present(coord)) then
       coord_usr = adjustl(coord)
       char1 = coord_usr(1:1)
       call ftupch(char1) ! uppercase
       if (char1 == 'C' .or. char1 == 'Q') then
          coord_usr = 'C'
          done_coord = .true.
       elseif (char1 == 'G') then
          coord_usr = 'G'
          done_coord = .true.
       elseif (char1 == 'E' ) then
          coord_usr='E'
          done_coord = .true.
       else
          print*,'Unrecognised COORD given : ',coord,' instead of C, G, or E'
          print*,'Proceed at you own risks '
          coord_usr = char1
          done_coord = .true.
       endif
       if (done_coord) then
          call ftukys(unit, 'COORDSYS',coord_usr,' Pixelisation coordinate system',status)
       endif
    endif


    !    write the user provided header literally, except for  PIXTYPE, TFORM*, TTYPE*, TUNIT*, INDXSCHM and GRAIN
    !    copy NSIDE, ORDERING and COORDSYS and POLAR if they are valid and not already given
    do i=1,nlheader
       card = header(i)
       if (card(1:5) == 'TTYPE' .or. card(1:5) == 'TFORM' .or. card(1:7) == 'PIXTYPE') then
          continue ! don't keep them
       else if (card(1:8) == 'INDXSCHM') then
          continue
       else if (card(1:5) == 'GRAIN') then ! already written above
          continue
       else if (card(1:13) == 'COMMENT GRAIN' .or. card(1:14) == 'COMMENT  GRAIN') then ! already written above
          continue
       else if (card(1:5) == 'TUNIT') then 
          if (trim(units_usr) == '') then
             call ftucrd(unit,'TUNIT'//card(6:7),card, status) !update TUNIT2 and above
          endif
       else if (card(1:5) == 'NSIDE') then
          call ftpsvc(card, svalue, comment, status)
          read(svalue,*) nside_hd
          npix_hd = nside2npix(nside_hd)
          if (.not. done_nside .and. npix_hd > 0) then
             call ftucrd(unit,'NSIDE',card, status) !update NSIDE
             done_nside = .true.
             nside_final = nside_hd
          endif
       else if (card(1:8) == 'ORDERING') then
          call ftpsvc(card, svalue, comment, status)
          svalue = adjustl(svalue)
          svalue = svalue(2:8) ! remove leading quote
          call ftupch(svalue)
          if (.not. done_order .and. (svalue(1:4)=='RING' .or. svalue(1:6) == 'NESTED')) then
             call ftucrd(unit,'ORDERING',card, status) !update ORDERING
             done_order = .true.
          endif
       else if (card(1:8) == 'COORDSYS') then
          if (.not. done_coord) call putrec(unit,card, status)
          done_coord = .true.
       else if (card(1:5) == 'POLAR') then
          if (.not. done_polar) then
             call ftucrd(unit,'POLAR', card, status) ! update POLAR
             done_polar = .true.
          endif
       else if (card(1:7) == 'EXTNAME') then
          call ftucrd(unit, 'EXTNAME', card, status)
       else if (card(1:1) /= ' ') then ! if none of the above, copy to FITS file
          call putrec(unit, card, status)
       endif
10     continue
    enddo


    ! test that required keywords have been provided in some way
    if (.not. done_nside) then
       print*,routine//'> NSIDE is a Healpix required keyword, '
       print*,routine//'>  it was NOT provided either as routine argument or in the input header'
       print*,routine//'>  abort execution, file not written'
       call fatal_error
    endif
    if (.not. done_order) then
       print*,routine//'> ORDER is a Healpix required keyword, '
       print*,routine//'>  it was NOT provided either as routine argument or in the input header'
       print*,routine//'>  abort execution, file not written'
       call fatal_error
    endif
!     if ((.not. done_polar) .and. extno_i >=2) then
!        print*,routine//'> Warning: POLAR keyword not set while 3 extensions have been written'
!     endif

    ! check that NSIDE is the same for all extensions
    if (extno_i == 0) then
       nside_old = nside_final
    else
       if (nside_final /= nside_old) then
          print*,routine//'> Inconsistent NSIDE: ',nside_final, nside_old
          print*,routine//'> Should use same NSIDE for all extensions'
       endif
    endif

    ! check validity of PIXEL
    npix_final = nside2npix(nside_final)
    if (kpix == I4B) then
       minpix = minval(pixel4(0:obs_npix8-1))
       maxpix = maxval(pixel4(0:obs_npix8-1))
    else
       minpix = minval(pixel8(0:obs_npix8-1))
       maxpix = maxval(pixel8(0:obs_npix8-1))
    endif
    if (minpix < 0 .or. maxpix > npix_final-1) then
       print*,routine//'> Actual pixel range = ',minpix,maxpix
       print*,routine//'> expected range (for Nside =',nside_final,') : 0, ',npix_final-1
       print*,routine//'> ABORT execution, file not written.'
       call fatal_error
    endif
    if (obs_npix8 > npix_final) then 
       print*,routine//'> The actual number of pixels ',obs_npix8
       print*,routine//'> is larger than the number of pixels over the whole sky : ',npix_final
       print*,routine//'> for Nside = ',nside_final
       print*,routine//'> ABORT execution, file not written.'
       call fatal_error
    endif

    !     write the extension one column by one column
    frow   = 1  ! starting position (row)
    felem  = 1  ! starting position (element)
    if (kpix == I4B) then
       call f90ftpcl_(unit, 1_i4b, frow, felem, obs_npix4, pixel4, status)
    else
       call f90ftpcl_(unit, 1_i4b, frow, felem, obs_npix4, pixel8, status)
    endif
    do i=1, nd2
       call f90ftpcl_(unit, i+1_i4b, frow, felem, obs_npix4, cutmap(0:,i), status)
    enddo

    !     ----------------------
    !     close and exit
    !     ----------------------

    !     close the file and free the unit number
    call ftclos(unit, status)

    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)

    return
  end subroutine sub_write_fits_partial_KLOAD
