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
module fitstools
  !
  ! module fitstools
  !  toward version 2.1
  !
  ! addition of read_fits_cut4 and write_fits_cut4
  ! improvement of input_map (addition of header output)
  ! improvement of read_asctab (more flexible on number of columns)
  ! removal of ask_* routines
  ! addition of get_clean_header
  ! addition of optional keyword units to input_map, read_fits_cut4, read_bintab
  ! 
  ! cleaning of the mapping between (l,m) <-> lm = l*l+l+m+1
  !  to remove round-off errors
  !
  ! v1.2_p01 : input_map, read_bintab : fix SUN compiler bug
  !            getsize_fits           : correct handling of ASCII table
  ! v1.2_p02   getsize_fits : correct handling of files with more than 2^31 elements
  !            read_bintab  : correction to array argument in ftgcve
  !            2004-Feb: 
  !               modification of dump_alms to deal with large arrays on MacOsX
  !               modification of npix in dump_alms
  !               modification of read_conbintab to deal with large arrrays
  !
  !  read_asctab : removal of nl_header from arguments list
  !
  ! 2004-07-21 : overload (SP/DP) output_map and write_bintab
  !   write_dbintab has been renamed write_plm
  ! 2004-11-01 : added function getnumext_fits
  !              added extno keyword to function getsize_fits
  ! 2004-12-22 : overload (SP/DP) read_asctab, dump_alms
  ! 2005-03-03 : added COORDSYS optional output to getsize_fits
  ! 2006-04-04 : edited to avoid concatenation problem in TFORMs (write_plm and write_bintabh with Ifort 9)
  ! 2007-02-16 : added EXTNO and POLARISATION keywords to Write_fits_cut4
  ! 2007-04-04 : close file on exit of getnumext_fits
  ! 2007-05-10 : increased maxdim (max number of columns per extension) from 20 to 40
  ! 2007-10-01 : getsize_fits returns -1 in case of error
  ! 2008-08-27 : in dump_alms and write_alms and write_*tab*: 
  !  do not write TTYPE# and TFORM# in excess of # of fields in the file
  ! 2010-11-23: implemented support for large (>2^31-1 pixel) map IO (ie, Nside > 8192)
  !             alm related IO is still limited to l<46340
  ! 2013-12-13 : increased MAXDIM from 40 to MAXDIM_TOP
  ! 2014-05-02: fixed problem with keywords having a long string value
  ! 2016-08-16: debugged input_map on cut4 FITS files with multi-HDU polarized data
  ! 2018-05-22: added unfold_weightsfile and unfold_weightslist
  ! -------------------------------------------------------------
  !
  ! --------------------------- from include file (see fits_template.f90)
  ! subroutine input_map
  ! subroutine read_bintab
  ! subroutine read_conbintab
  ! subroutine write_bintab (write_bintab4 and write_bintab8)
  ! subroutine write_asctab
  ! subroutine dump_alms
  ! subroutine write_alms
  ! subroutine read_alms
  ! subroutine read_bintod
  ! subroutine write_bintabh
  ! subroutine unfold_weights
  ! subroutine read_fits_partial
  ! subroutine write_fits_partial
  ! ----------------------------------
  !
  ! subroutine read_fits_cut4               ?
  ! subroutine write_fits_cut4              ?
  ! subroutine fits2cl                       NA    accepts extname addressing (eg myfile.fits[beam])
  ! subroutine read_asctab                   NA
  !
  ! subroutine output_map                  (4/8)
  ! subroutine write_dbintab   : Obsolete
  ! subroutine write_plm                    NA
  !
  ! subroutine alms2fits                    NotYet
  ! subroutine fits2alms                    NotYet
  ! subroutine input_tod                   (8)
  ! subroutine read_dbintab                NA
  !                                       
  ! subroutine printerror                  NA
  !                                       
  ! subroutine read_par                    NA
  ! function getnumext_fits                NA
  ! function getsize_fits                  (8)    accepts extname addressing (eg myfile.fits[beam])
  ! function number_of_alms                NotYet
  !                                       
  ! subroutine putrec                      NA
  ! subroutine get_clean_header            NA
  !
  ! subroutine check_input_map
  !
  use healpix_types
  use misc_utils, only : fatal_error, assert, assert_alloc, strupcase, string
  implicit none

  real(kind=SP),     private, parameter :: s_bad_value = HPX_SBADVAL
  real(kind=DP),     private, parameter :: d_bad_value = HPX_DBADVAL
  integer(kind=I4B), private, parameter :: i_bad_value = -1637500000
  integer(kind=I8B), private, parameter :: l_bad_value = -1637500000_I8B
  integer(I4B) ,     private, parameter :: nchunk_max  = 12000
  integer(I4B),      private, parameter :: MAXDIM_TOP  = 199 ! < 999

!   interface read_fits_cut4
!      module procedure read_fits_cut4_s,read_fits_cut4_d
!   end interface

!   interface write_fits_cut4
!      module procedure write_fits_cut4_s,write_fits_cut4_d
!   end interface


  interface input_map
#ifdef NO64BITS
     module procedure input_map8_s,input_map8_d
#else
     module procedure input_map8_s,input_map8_d, input_map4_s,input_map4_d
#endif
  end interface

  interface output_map
     module procedure output_map_s, output_map_d
  end interface

  interface write_plm
#ifdef NO64BITS
     module procedure write_plm8
#else
     module procedure write_plm8, write_plm4
#endif
  end interface

  interface write_bintab
#ifdef NO64BITS
     module procedure write_bintab8_s, write_bintab8_d
#else
     module procedure write_bintab8_s, write_bintab8_d, write_bintab4_s, write_bintab4_d
#endif
  end interface

  interface read_bintab
#ifdef NO64BITS
     module procedure read_bintab8_s, read_bintab8_d
#else
     module procedure read_bintab8_s, read_bintab8_d, read_bintab4_s, read_bintab4_d
#endif
  end interface

  interface read_fits_partial
#ifdef NO64BITS
     module procedure read_fits_partial8_s, read_fits_partial8_d
#else
     module procedure read_fits_partial8_s, read_fits_partial8_d, &
          & read_fits_partial4_s, read_fits_partial4_d
#endif
  end interface

  interface write_fits_partial
#ifdef NO64BITS
     module procedure write_fits_partial4_s, write_fits_partial4_d
#else
     module procedure write_fits_partial8_s, write_fits_partial8_d, &
          & write_fits_partial4_s, write_fits_partial4_d
#endif
  end interface

  !
  interface input_tod
     module procedure input_tod_s,input_tod_d
  end interface

  interface read_bintod
     module procedure read_bintod_s,read_bintod_d
  end interface

  interface write_bintabh
     module procedure write_bintabh_s,write_bintabh_d
  end interface
  !
  interface write_asctab
     module procedure write_asctab_s, write_asctab_d
  end interface

  interface fits2cl
     module procedure fits2cl_s, fits2cl_d
  end interface

  interface read_asctab
     module procedure read_asctab_s, read_asctab_d, read_asctab_v12s, read_asctab_v12d
  end interface
  !
  interface fits2alms
     module procedure fits2alms_s, fits2alms_d
  end interface

  interface alms2fits
     module procedure alms2fits_s, alms2fits_d
  end interface

  interface read_conbintab
     module procedure read_conbintab_s, read_conbintab_d
  end interface

  interface dump_alms
     module procedure dump_alms_s, dump_alms_d
  end interface

  interface read_alms
     module procedure read_alms_s, read_alms_d
  end interface

  interface write_alms
     module procedure write_alms_s, write_alms_d
  end interface

  interface f90ftpcl_
#ifdef NO64BITS
     module procedure f90ftpcle, f90ftpcld, f90ftpclj
#else
     module procedure f90ftpcle, f90ftpcld, f90ftpclj, f90ftpclk
#endif
  end interface

  interface f90ftgcv_
#ifdef NO64BITS
     module procedure f90ftgcve, f90ftgcvd, f90ftgcvj
#else
     module procedure f90ftgcve, f90ftgcvd, f90ftgcvj, f90ftgcvk
#endif
  end interface

  interface f90ftgpv_
     module procedure f90ftgpve, f90ftgpvd
  end interface

  interface f90ftgky_
#ifdef NO64BITS
     module procedure f90ftgkye, f90ftgkyd, f90ftgkyj
#else
     module procedure f90ftgkye, f90ftgkyd, f90ftgkyj, f90ftgkyk
#endif
  end interface

  interface map_bad_pixels
     module procedure map_bad_pixels_s, map_bad_pixels_d
  end interface

  interface unfold_weightsfile
     module procedure unfold_weightsfile_s, unfold_weightsfile_d
  end interface

  interface unfold_weightslist
     module procedure unfold_weightslist_s, unfold_weightslist_d
  end interface

  private

  public :: read_fits_cut4, write_fits_cut4, & 
       & input_map, read_bintab,  &
       & output_map, write_bintab
  public :: read_fits_partial, write_fits_partial
  public :: write_plm
  public :: fits2cl, read_asctab, write_asctab
  public :: read_dbintab, write_dbintab
  public :: input_tod, write_bintabh
  public :: dump_alms, alms2fits, fits2alms, write_alms, read_alms, read_conbintab
  public :: printerror, read_par, getsize_fits, number_of_alms, getnumext_fits
  public :: putrec
  public :: check_input_map
  public :: unfold_weightsfile, unfold_weightslist

contains

  !-------------------------------------------------------------------------------
  ! generic interface F90FTPCL_ for FITSIO's FTPCL[E,D,J,K]
  !           writes data in ASCTAB or BINTAB
  subroutine f90ftpcle(unit, colnum, frow, felem, np, data, status)
    integer(I4B), intent(in)  :: unit, colnum, frow, felem, np
    integer(I4B), intent(out) :: status
    real(SP),     intent(in), dimension(0:)  :: data
    call ftpcle(unit, colnum, frow, felem, np, data, status)
    return
  end subroutine f90ftpcle
  subroutine f90ftpcld(unit, colnum, frow, felem, np, data, status)
    integer(I4B), intent(in)  :: unit, colnum, frow, felem, np
    integer(I4B), intent(out) :: status
    real(DP),     intent(in), dimension(0:)  :: data
    call ftpcld(unit, colnum, frow, felem, np, data, status)
    return
  end subroutine f90ftpcld
  subroutine f90ftpclj(unit, colnum, frow, felem, np, data, status)
    integer(I4B), intent(in)  :: unit, colnum, frow, felem, np
    integer(I4B), intent(out) :: status
    integer(I4B), intent(in), dimension(0:)  :: data
    call ftpclj(unit, colnum, frow, felem, np, data, status)
    return
  end subroutine f90ftpclj
#ifndef NO64BITS
  subroutine f90ftpclk(unit, colnum, frow, felem, np, data, status)
    integer(I4B), intent(in)  :: unit, colnum, frow, felem, np
    integer(I4B), intent(out) :: status
    integer(I8B), intent(in), dimension(0:)  :: data
    call ftpclk(unit, colnum, frow, felem, np, data, status)
    return
  end subroutine f90ftpclk
#endif
  !-------------------------------------------------------------------------------
  ! generic interface F90FTGCV_ for FITSIO's FTGCV[E,D,J,K]
  !           reads data from BINTAB
  subroutine f90ftgcve(unit, colnum, frow, felem, np, nullval, data, anynull, status)
    integer(I4B), intent(in)  :: unit, colnum, frow, felem, np
    integer(I4B), intent(out) :: status
    logical(LGT), intent(out) :: anynull
    real(SP),     intent(out), dimension(0:) :: data
    real(SP),     intent(in)                 :: nullval
    call ftgcve(unit, colnum, frow, felem, np, nullval, data, anynull, status)
    return
  end subroutine f90ftgcve
  subroutine f90ftgcvd(unit, colnum, frow, felem, np, nullval, data, anynull, status)
    integer(I4B), intent(in)  :: unit, colnum, frow, felem, np
    integer(I4B), intent(out) :: status
    logical(LGT), intent(out) :: anynull
    real(DP),     intent(out), dimension(0:) :: data
    real(DP),     intent(in)                 :: nullval
    call ftgcvd(unit, colnum, frow, felem, np, nullval, data, anynull, status)
    return
  end subroutine f90ftgcvd
  subroutine f90ftgcvj(unit, colnum, frow, felem, np, nullval, data, anynull, status)
    integer(I4B), intent(in)  :: unit, colnum, frow, felem, np
    integer(I4B), intent(out) :: status
    logical(LGT), intent(out) :: anynull
    integer(I4B), intent(out), dimension(0:) :: data
    integer(I4B), intent(in)                 :: nullval
    call ftgcvj(unit, colnum, frow, felem, np, nullval, data, anynull, status)
    return
  end subroutine f90ftgcvj
#ifndef NO64BITS
  subroutine f90ftgcvk(unit, colnum, frow, felem, np, nullval, data, anynull, status)
    integer(I4B), intent(in)  :: unit, colnum, frow, felem, np
    integer(I4B), intent(out) :: status
    logical(LGT), intent(out) :: anynull
    integer(I8B), intent(out), dimension(0:) :: data
    integer(I8B), intent(in)                 :: nullval
    call ftgcvk(unit, colnum, frow, felem, np, nullval, data, anynull, status)
    return
  end subroutine f90ftgcvk
#endif
  !-------------------------------------------------------------------------------
  ! generic interface F90FTGPV_ for FITSIO's FTGPVE and FTGPVD
  !           reads data from IMAGE
  subroutine f90ftgpve(unit, group, firstpix, np, nullval, data, anynull, status)
    integer(I4B), intent(in)  :: unit, group, firstpix, np
    integer(I4B), intent(out) :: status
    logical(LGT), intent(out) :: anynull
    real(SP),     intent(out), dimension(0:) :: data
    real(SP),     intent(in)                 :: nullval
    call ftgpve(unit, group, firstpix, np, nullval, data, anynull, status)
    return
  end subroutine f90ftgpve
  subroutine f90ftgpvd(unit, group, firstpix, np, nullval, data, anynull, status)
    integer(I4B), intent(in)  :: unit, group, firstpix, np
    integer(I4B), intent(out) :: status
    logical(LGT), intent(out) :: anynull
    real(DP),     intent(out), dimension(0:) :: data
    real(DP),     intent(in)                 :: nullval
    call ftgpvd(unit, group, firstpix, np, nullval, data, anynull, status)
    return
  end subroutine f90ftgpvd
  !-------------------------------------------------------------------------------
  ! generic interface F90FTGKY_ for FITSIO's FTGKY[E,D,J,K]
  !           reads a keyword
  subroutine f90ftgkye(unit, keyword, value, comment, status)
    integer(I4B),     intent(in)  :: unit
    character(len=*), intent(in)  :: keyword
    integer(I4B),     intent(out) :: status
    character(len=*), intent(out) :: comment
    real(SP),         intent(out) :: value
    call ftgkye(unit, keyword, value, comment, status)
    return
  end subroutine f90ftgkye
  subroutine f90ftgkyd(unit, keyword, value, comment, status)
    integer(I4B),     intent(in)  :: unit
    character(len=*), intent(in)  :: keyword
    integer(I4B),     intent(out) :: status
    character(len=*), intent(out) :: comment
    real(DP),         intent(out) :: value
    call ftgkyd(unit, keyword, value, comment, status)
    return
  end subroutine f90ftgkyd
  subroutine f90ftgkyj(unit, keyword, value, comment, status)
    integer(I4B),     intent(in)  :: unit
    character(len=*), intent(in)  :: keyword
    integer(I4B),     intent(out) :: status
    character(len=*), intent(out) :: comment
    integer(I4B),     intent(out) :: value
    call ftgkyj(unit, keyword, value, comment, status)
    return
  end subroutine f90ftgkyj
#ifndef NO64BITS
  subroutine f90ftgkyk(unit, keyword, value, comment, status)
    integer(I4B),     intent(in)  :: unit
    character(len=*), intent(in)  :: keyword
    integer(I4B),     intent(out) :: status
    character(len=*), intent(out) :: comment
    integer(I8B),     intent(out) :: value
    call ftgkyk(unit, keyword, value, comment, status)
    return
  end subroutine f90ftgkyk
#endif
  !-------------------------------------------------------------------------------


  ! define routine with SP I/O
#include "fits_s_inc.F90"

  ! define routine with DP I/O
#include "fits_d_inc.F90"



  !=======================================================================
  subroutine read_fits_cut4(filename, npixtot, pixel, signal, n_obs, serror, header, units, extno)
    !=======================================================================
    !
    ! routine to read FITS file with 4 columns : PIXEL, SIGNAL, N_OBS, SERROR
    !
    ! read_fits_cut4(filename, npixtot, pixel, 
    !      [signal, n_obs, serror, header, units, extno])
    !=======================================================================
    integer(I4b), parameter :: KMAP = SP

    character(len=*),                   intent(in)              :: filename
    integer(I4B),                       intent(in)              :: npixtot
    integer(I4B), dimension(0:npixtot-1), intent(out)           :: pixel
    real(KMAP),   dimension(0:npixtot-1), intent(out), optional :: signal
    integer(I4B), dimension(0:npixtot-1), intent(out), optional :: n_obs
    real(KMAP),   dimension(0:npixtot-1), intent(out), optional :: serror
    character(len=*), dimension(1:),      intent(out), optional :: header
    character(len=*),                     intent(out), optional :: units
    integer(I4B),                         intent(in),  optional :: extno

    integer(I4B) :: obs_npix

    integer(I4B), parameter :: MAXDIM = MAXDIM_TOP !number of columns in the extension
    integer(I4B) :: blocksize, datacode
    integer(I4B) :: firstpix, frow, hdutype
    integer(I4B) :: naxis, nfound, nmove, npix, nrows
    integer(I4B) :: readwrite, status, tfields, unit, varidat, width
    integer(I4B) :: repeat, repeat1, repeat2
    logical(LGT) :: anynull, extend
    character(len=20), dimension(MAXDIM) :: ttype, tform, tunit
    character(len=20) :: extname
    character(len=80) :: comment
    real(KMAP) :: nullval
    !=======================================================================

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
       print*,'an image was found in the FITS file '//trim(filename)
       print*,'... it is ignored.'
    endif

    !     determines the presence of an extension
    call ftgkyl(unit,'EXTEND', extend, comment, status)
    if (status > 0) then 
       print*,'extension expected and not found in FITS file '//trim(filename)
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
    if (tfields < 4) then
       print*,'Expected 4 columns in FITS file '//trim(filename)
       print*,'found ',tfields
       if (tfields < 2) call fatal_error
       if (.not.  (trim(ttype(1)) == 'PIXEL' &
            & .or. trim(ttype(1)) == 'PIX'  ) ) call fatal_error
    endif

    if (present(header)) then 
       header = ""
       status = 0
       call get_clean_header(unit, header, filename, status)
    endif

    if (present(units)) then 
       ! the second column contains the SIGNAL, for which we 
       ! already read the units from the header
       units = adjustl(tunit(2))
    endif

    !        finds the bad data value
!     if (KMAP == SP) CALL ftgkye(unit,'BAD_DATA',nullval,comment,status)
!     if (KMAP == DP) CALL ftgkyd(unit,'BAD_DATA',nullval,comment,status)
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

    call ftgkyj(unit,'OBS_NPIX',obs_npix,comment,status)
    if (status == 202) then ! obs_npix not found
       obs_npix = nrows * repeat
       status = 0
    endif

    npix = min(npixtot, obs_npix)
    call ftgcvj(unit, 1_i4b, frow, firstpix, npix, i_bad_value, pixel(0), anynull, status)
    if (present(signal)) then
       call f90ftgcv_(unit, 2_i4b, frow, firstpix, npix, nullval, signal(0:npix-1), anynull, status)
    endif
    if (tfields >= 3 .and. present(n_obs)) &
         & call ftgcvj(unit, 3_i4b, frow, firstpix, npix, i_bad_value, n_obs(0), anynull, status)
    if (tfields >= 4 .and. present(serror)) then
       call f90ftgcv_(unit, 4_i4b, frow, firstpix, npix, nullval, serror(0:npix-1), anynull, status)
    endif

    !     close the file
    call ftclos(unit, status)

    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)

    return
  end subroutine read_fits_cut4

  !=======================================================================
  subroutine write_fits_cut4(filename, obs_npix, pixel, signal, n_obs, serror, &
       &                     header, coord, nside, order, units, extno, polarisation)
    !=======================================================================
    ! writes a fits file for cut sky data set with 4 columns:
    !  PIXEL, SIGNAL, N_OBS, SERROR
    ! 
    ! write_fits_cut4(filename, obs_npix, pixel, signal, n_obs, serror
    !         [, header, coord, nside, order, units, extno, polarisation])
    !=======================================================================
    use pix_tools, only: nside2npix
    integer(I4b), parameter :: KMAP = SP

    character(len=*), intent(in) :: filename
    integer(I4B)     :: obs_npix
    integer(I4B),  dimension(0:obs_npix-1), intent(in)           :: pixel
    real(KMAP),    dimension(0:obs_npix-1), intent(in)           :: signal
    integer(I4B),  dimension(0:obs_npix-1), intent(in)           :: n_obs
    real(KMAP),    dimension(0:obs_npix-1), intent(in)           :: serror
    character(len=*), dimension(1:),       intent(in), optional :: header
    integer(I4B),                          intent(in), optional :: nside, order
    character(len=*),                      intent(in), optional :: coord, units
    integer(I4B),                          intent(in), optional :: extno, polarisation

    character(len=*), parameter :: routine = 'write_fits_cut4'
    ! --- healpix/fits related variables ----
    integer(I4B)     :: ncol, grain
    integer(I4B)     :: npix_hd, nside_hd, nside_final, npix_final, i
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
    character(len=1) :: pform
    character(len=filenamelen) sfilename

    integer(I4B), save :: nside_old
    !=======================================================================
    if (KMAP == SP) pform='E'
    if (KMAP == DP) pform='D'

    extnames = (/ 'TEMPERATURE   ' , & ! default in case of polarisation
         &        'Q_POLARISATION' , &
         &        'U_POLARISATION' /)

    ncol = 4
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
    if (present(polarisation)) then
       polar_flag = (polarisation /= 0)
       done_polar = .true. ! will ignore the header provided value of POLAR
    endif

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
    if (obs_npix < repeat) repeat = 1
    nrows    = (obs_npix + repeat - 1)/ repeat ! naxis1
    write(srepeat,'(i4)') repeat
    srepeat = adjustl(srepeat)

    repeatg = repeat * grain
    write(srepeatg,'(i4)') repeatg
    srepeatg = adjustl(srepeatg)

    tfields  = ncol
    ttype(1:4) = (/ 'PIXEL ','SIGNAL','N_OBS ','SERROR' /)
    tform(1) = trim(srepeat)//'J'
    tform(2) = trim(srepeatg)//pform
    tform(3) = trim(srepeatg)//'J'
    tform(4) = trim(srepeatg)//pform

    tunit =  ' '      ! optional, will not appear
    tunit(2) = units_usr
    tunit(4) = units_usr
    extname  = 'SKY_OBSERVATION'      ! default, will be overridden by user provided one if any
    if (polar_flag) extname = extnames(1+extno_i)
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
    call ftpkyj(unit,'OBS_NPIX',obs_npix, ' Number of pixels observed and recorded',status)


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
          if ((card(6:6) == '2' .or. card(6:6) == '4') .and. trim(units_usr) == '') then
             call ftucrd(unit,'TUNIT'//card(6:6),card, status) !update TUNIT2 and TUNIT4
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
    if ((.not. done_polar) .and. extno_i >=2) then
       print*,routine//'> Warning: POLAR keyword not set while 3 extensions have been written'
    endif

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
    minpix = minval(pixel(0:obs_npix-1))
    maxpix = maxval(pixel(0:obs_npix-1))
    if (minpix < 0 .or. maxpix > npix_final-1) then
       print*,routine//'> Actual pixel range = ',minpix,maxpix
       print*,routine//'> expected range (for Nside =',nside_final,') : 0, ',npix_final-1
       print*,routine//'> ABORT execution, file not written.'
       call fatal_error
    endif
    if (obs_npix > npix_final) then 
       print*,routine//'> The actual number of pixels ',obs_npix
       print*,routine//'> is larger than the number of pixels over the whole sky : ',npix_final
       print*,routine//'> for Nside = ',nside_final
       print*,routine//'> ABORT execution, file not written.'
       call fatal_error
    endif

    !     write the extension one column by one column
    frow   = 1  ! starting position (row)
    felem  = 1  ! starting position (element)
    call ftpclj(unit, 1_i4b, frow, felem, obs_npix, pixel , status)
    call f90ftpcl_(unit, 2_i4b, frow, felem, obs_npix, signal, status)

    if (tfields >= 3) call ftpclj(unit, 3_i4b, frow, felem, obs_npix, n_obs , status)
    if (tfields >= 4) then
       call f90ftpcl_(unit, 4_i4b, frow, felem, obs_npix, serror, status)
    endif

    !     ----------------------
    !     close and exit
    !     ----------------------

    !     close the file and free the unit number
    call ftclos(unit, status)

    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)

    return
  end subroutine write_fits_cut4

  !=======================================================================
  ! FITS2CL(filename, clin, lmax, ncl, header, units, fmissval)
  !     Read C_l from a FITS file
  !   Aug 2000 : modification by EH of the number of columns actually read
  !
  !     This routine is used for reading input power spectra for synfast
  !
  !  Dec 2004: overloading for single and double precision output array (clin)
  ! Feb 2013: added fmissval (value to be given to NaN or Inf C(l))
  !      if absent, those C(l) are left unchanged
  !   
  !=======================================================================
  subroutine fits2cl_s(filename, clin, lmax, ncl, header, units, fmissval)
    !=======================================================================
    !        single precision
    !=======================================================================
    integer(I4b), parameter :: KCL = SP
    !
    CHARACTER(LEN=*),                          INTENT(IN) :: filename
    INTEGER(I4B),                              INTENT(IN) :: lmax, ncl
    REAL(KCL),        DIMENSION(0:lmax,1:ncl), INTENT(OUT) :: clin
    CHARACTER(LEN=*), DIMENSION(1:),            INTENT(OUT) :: header
    CHARACTER(LEN=*), dimension(1:), optional,  INTENT(OUT) :: units
    real(KCL),                       optional,  intent(IN):: fmissval

    real(DP), dimension(:,:), allocatable :: cl_dp
    real(DP) :: fmvd

    ! since the arrays involved are small
    ! read in double precision, and then cast in single
    allocate(cl_dp(0:lmax, 1:ncl))
    if (present(fmissval)) then
       fmvd = fmissval * 1.0_DP
       call fits2cl_d(filename, cl_dp, lmax, ncl, header, units, fmvd)
    else
       call fits2cl_d(filename, cl_dp, lmax, ncl, header, units)
    endif
    clin(0:lmax, 1:ncl) = cl_dp(0:lmax, 1:ncl)

    return
  end subroutine fits2cl_s

  subroutine fits2cl_d(filename, clin, lmax, ncl, header, units, fmissval)
    !=======================================================================
    !        double precision
    !=======================================================================
    integer(I4b), parameter :: KCL = DP
    !
    CHARACTER(LEN=*),                          INTENT(IN) :: filename
    INTEGER(I4B),                              INTENT(IN) :: lmax, ncl
    REAL(KCL),         DIMENSION(0:lmax, 1:ncl), INTENT(OUT) :: clin
    CHARACTER(LEN=*), DIMENSION(1:),   INTENT(OUT) :: header
    CHARACTER(LEN=*), dimension(1:), optional,  INTENT(OUT) :: units
    real(KCL),                       optional,  intent(IN):: fmissval

    INTEGER(I4B) :: status,unit,readwrite,blocksize,naxes(2),ncl_file, naxis
    INTEGER(I4B) :: firstpix,lmax_file,lmax_min,nelems,datacode,repeat,width
    CHARACTER(LEN=80) :: comment, pdmstr
    LOGICAL :: extend
    INTEGER(I4B) :: nmove, hdutype, hdunum
    INTEGER(I4B) :: column, frow
    REAL(KCL), DIMENSION(:), ALLOCATABLE :: clin_file

    INTEGER(I4B), PARAMETER :: MAXDIM = MAXDIM_TOP !number of columns in the extension
    INTEGER(I4B) :: rowlen, nrows, varidat
    INTEGER(I4B),      dimension(1:MAXDIM) :: tbcol
    CHARACTER(LEN=20), dimension(1:MAXDIM) :: ttype, tform, tunit
    CHARACTER(LEN=20)                      :: extname
    LOGICAL :: anynull
    REAL(KCL) ::  nullval
    logical :: planck_format
    integer(i8b) :: nbads


    !-----------------------------------------------------------------------
    status=0
    anynull = .false.
    ttype=''
    tform=''
    tunit=''
    extname=''
    comment=''

    unit = 110
    naxes(1) = 1
    naxes(2) = 1

    readwrite=0
    call ftnopn(unit,filename,readwrite, status) !open primary HDU or specified HDU
    if (status > 0) call printerror(status)
    !     -----------------------------------------
    call ftghdn(unit, hdunum)
    if (hdunum == 1) then  ! in primary HDU: move to next HDU
       !     determines the presence of image
       call ftgkyj(unit,'NAXIS', naxis, comment, status)
       if (status > 0) call printerror(status)

       !     determines the presence of an extension
       call ftgkyl(unit,'EXTEND', extend, comment, status)
       if (status > 0) call printerror(status)

       nmove = +1
       call ftmrhd(unit, nmove, hdutype, status)
    else ! already in non primary HDU
       extend = .true.
       call ftghdt(unit, hdutype, status)
    endif

    if (extend) then ! there is an extension

       call assert ((hdutype==1).or.(hdutype==2), 'this is not a table')

       !        reads keywords related to table layout
       if (hdutype==1) then ! ASCII table
         call ftghtb(unit, MAXDIM, rowlen, &
              &      nrows, ncl_file, ttype, tbcol, &
              &      tform, tunit, extname, status)
         repeat = 1
       else ! binary table
         call ftghbn(unit, MAXDIM, &
              &      nrows, ncl_file, ttype,        &
              &      tform, tunit,extname, varidat, status)
          call ftbnfm(tform(1), datacode, repeat, width, status)
       endif

       status = 0
       header = ""
       call get_clean_header(unit, header, filename, status)

       if (present(units)) then 
          do column = 1, min(ncl, ncl_file)
             units(column) = adjustl(tunit(column))
          enddo
       endif

       !        reads the columns
       column = 1
       frow = 1
       firstpix = 1
       lmax_file = nrows*repeat - 1
       lmax_min = MIN(lmax,lmax_file)
       nelems = lmax_min + 1
       nullval = 0.0_KCL ! CFITSIO will leave bad pixels unchanged, and so will this routine
       if (present(fmissval)) then
          if (fmissval == 0.0_KCL) then
             nullval = HPX_DBADVAL ! sentinel value to be given to bad pixels by CFITSIO, will later be mapped to user defined 0
          else
             nullval = fmissval ! user defined value to be given to bad pixels by CFITSIO, and by this routine
          endif
       endif

! check for the special Planck format (i.e. one additional column)
       planck_format=.true.
       call ftgkys(unit,'PDMTYPE',pdmstr,comment,status)
       if (status==202) then
         planck_format=.false.
         status=0
       endif
       allocate(clin_file(0:lmax_min),stat=status)
       clin = 0.0_KCL                         ! modification by EH
       if (planck_format) then
         do column = 1, MIN(ncl, ncl_file-1) ! modification by EH
            clin_file(:) = 0.0_KCL
            call ftgcvd(unit, column+1_i4b, frow, firstpix, nelems, nullval, &
                 &        clin_file(0), anynull, status)
            clin(0:lmax_min, column) = clin_file(0:lmax_min)
         enddo
       else
         do column = 1, MIN(ncl, ncl_file) ! modification by EH
            clin_file(:) = 0.0_KCL
            call ftgcvd(unit, column, frow, firstpix, nelems, nullval, &
                 &        clin_file(0), anynull, status)
            clin(0:lmax_min, column) = clin_file(0:lmax_min)
         enddo
       endif
       deallocate(clin_file)
       if (present(fmissval)) then
          if (nullval /= fmissval) call map_bad_pixels(clin, nullval, fmissval, nbads)
       endif
    else ! no image no extension, you are dead, man
       call fatal_error(' No image, no extension in '//trim(filename))
    endif

    !     close the file
    call ftclos(unit, status)

    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)

!     CHARACTER(LEN=*),                          INTENT(IN) :: filename
!     INTEGER(I4B),                              INTENT(IN) :: lmax, ncl
!     REAL(KCL),         DIMENSION(0:lmax, 1:ncl), INTENT(OUT) :: clin
!     CHARACTER(LEN=*), DIMENSION(1:),   INTENT(OUT) :: header
!     CHARACTER(LEN=*), dimension(1:), optional,  INTENT(OUT) :: units

!     INTEGER(I4B) :: status,unit,readwrite,blocksize,naxes(2),ncl_file, naxis
!     INTEGER(I4B) :: firstpix,lmax_file,lmax_min,nelems
!     CHARACTER(LEN=80) :: comment, pdmstr
!     LOGICAL :: extend
!     INTEGER(I4B) :: nmove, hdutype
!     INTEGER(I4B) :: column, frow
!     REAL(KCL), DIMENSION(:), ALLOCATABLE :: clin_file

!     INTEGER(I4B), PARAMETER :: MAXDIM = MAXDIM_TOP !number of columns in the extension
!     INTEGER(I4B) :: rowlen, nrows, varidat
!     INTEGER(I4B),      dimension(1:MAXDIM) :: tbcol
!     CHARACTER(LEN=20), dimension(1:MAXDIM) :: ttype, tform, tunit
!     CHARACTER(LEN=20)                      :: extname
!     LOGICAL :: anynull
!     REAL(KCL) ::  nullval
!     logical :: planck_format


!     !-----------------------------------------------------------------------
!     status=0
!     anynull = .false.
!     ttype=''
!     tform=''
!     tunit=''
!     extname=''
!     comment=''

!     unit = 110
!     naxes(1) = 1
!     naxes(2) = 1

!     readwrite=0
!     call ftopen(unit,filename,readwrite,blocksize,status)
!     if (status > 0) call printerror(status)
!     !     -----------------------------------------

!     !     determines the presence of image
!     call ftgkyj(unit,'NAXIS', naxis, comment, status)
!     if (status > 0) call printerror(status)

!     !     determines the presence of an extension
!     call ftgkyl(unit,'EXTEND', extend, comment, status)
!     if (status > 0) call printerror(status)

!     if (extend) then ! there is an extension
!        nmove = +1
!        call ftmrhd(unit, nmove, hdutype, status)

!        call assert ((hdutype==1).or.(hdutype==2), 'this is not a table')

!        !        reads keywords related to table layout
!        if (hdutype==1) then
!          call ftghtb(unit, MAXDIM, &
!               &        rowlen, nrows, ncl_file, ttype, tbcol, tform, tunit, &
!               &        extname, status)
!        else
!          call ftghbn(unit, MAXDIM, &
!               &        nrows, ncl_file, ttype, tform, tunit,extname, &
!               &        varidat, status)
!        endif

!        status = 0
!        header = ""
!        call get_clean_header(unit, header, filename, status)

!        if (present(units)) then 
!           do column = 1, min(ncl, ncl_file)
!              units(column) = adjustl(tunit(column))
!           enddo
!        endif

!        !        reads the columns
!        column = 1
!        frow = 1
!        firstpix = 1
!        lmax_file = nrows - 1
!        lmax_min = MIN(lmax,lmax_file)
!        nullval = 0.0_KCL
!        nelems = lmax_min + 1

! ! check for the special Planck format (i.e. one additional column)
!        planck_format=.true.
!        call ftgkys(unit,'PDMTYPE',pdmstr,comment,status)
!        if (status==202) then
!          planck_format=.false.
!          status=0
!        endif
!        allocate(clin_file(0:lmax_min),stat=status)
!        clin = 0.0_KCL                         ! modification by EH
!        if (planck_format) then
!          do column = 1, MIN(ncl, ncl_file-1) ! modification by EH
!             clin_file(:) = 0.0_KCL
!             call ftgcvd(unit, column+1_i4b, frow, firstpix, nelems, nullval, &
!                  &        clin_file(0), anynull, status)
!             clin(0:lmax_min, column) = clin_file(0:lmax_min)
!          enddo
!        else
!          do column = 1, MIN(ncl, ncl_file) ! modification by EH
!             clin_file(:) = 0.0_KCL
!             call ftgcvd(unit, column, frow, firstpix, nelems, nullval, &
!                  &        clin_file(0), anynull, status)
!             clin(0:lmax_min, column) = clin_file(0:lmax_min)
!          enddo
!        endif
!        deallocate(clin_file)
!     else ! no image no extension, you are dead, man
!        call fatal_error(' No image, no extension')
!     endif

!     !     close the file
!     call ftclos(unit, status)

!     !     check for any error, and if so print out error messages
!     if (status > 0) call printerror(status)



    return
  end subroutine fits2cl_d

  !=======================================================================
  ! READ_ASC_TAB(filename, clin, lmax, ncl, header, units)
  ! superseded by fits2cl
  !=======================================================================
  subroutine read_asctab_v12s(filename, clin, lmax, ncl, header, nlheader, units)
    !=======================================================================
    !        single precision, legacy code
    !=======================================================================
    integer(I4b), parameter :: KCL = SP
    !
    CHARACTER(LEN=*),                          INTENT(IN) :: filename
    INTEGER(I4B),                              INTENT(IN) :: lmax, ncl,nlheader
    REAL(KCL),         DIMENSION(0:lmax,1:ncl), INTENT(OUT) :: clin
    CHARACTER(LEN=*), DIMENSION(1:),   INTENT(OUT) :: header
    CHARACTER(LEN=*), dimension(1:), optional,  INTENT(OUT) :: units

    print*,"-------------------------------------------------------------"
    print*,"WARNING : the routine read_asctab is now obsolete"
    print*,"  Use "
    print*," call fits2cl(filename, clin, lmax, ncl, header, units)"
    print*,"  instead (same module)"
    print*,"-------------------------------------------------------------"

    call fits2cl(filename, clin, lmax, ncl, header, units)
    return
  end subroutine read_asctab_v12s
  subroutine read_asctab_v12d(filename, clin, lmax, ncl, header, nlheader, units)
    !=======================================================================
    !        double precision, legacy code
    !=======================================================================
    integer(I4b), parameter :: KCL = DP
    !
    CHARACTER(LEN=*),                          INTENT(IN) :: filename
    INTEGER(I4B),                              INTENT(IN) :: lmax, ncl,nlheader
    REAL(KCL),         DIMENSION(0:lmax,1:ncl), INTENT(OUT) :: clin
    CHARACTER(LEN=*), DIMENSION(1:),   INTENT(OUT) :: header
    CHARACTER(LEN=*), dimension(1:), optional,  INTENT(OUT) :: units

    print*,"-------------------------------------------------------------"
    print*,"WARNING : the routine read_asctab is now obsolete"
    print*,"  Use "
    print*," call fits2cl(filename, clin, lmax, ncl, header, units)"
    print*,"  instead (same module)"
    print*,"-------------------------------------------------------------"

    call fits2cl(filename, clin, lmax, ncl, header, units)
    return
  end subroutine read_asctab_v12d
  !=======================================================================
  subroutine read_asctab_s(filename, clin, lmax, ncl, header, units)
    !=======================================================================
    !        single precision
    !=======================================================================
    integer(I4b), parameter :: KCL = SP
    !
    CHARACTER(LEN=*),                          INTENT(IN) :: filename
    INTEGER(I4B),                              INTENT(IN) :: lmax, ncl
    REAL(KCL),         DIMENSION(0:lmax,1:ncl), INTENT(OUT) :: clin
    CHARACTER(LEN=*), DIMENSION(1:),   INTENT(OUT) :: header
    CHARACTER(LEN=*), dimension(1:), optional,  INTENT(OUT) :: units

    print*,"-------------------------------------------------------------"
    print*,"WARNING : the routine read_asctab is now obsolete"
    print*,"  Use "
    print*," call fits2cl(filename, clin, lmax, ncl, header, units)"
    print*,"  instead (same module)"
    print*,"-------------------------------------------------------------"

    call fits2cl(filename, clin, lmax, ncl, header, units)
    return
  end subroutine read_asctab_s

  subroutine read_asctab_d(filename, clin, lmax, ncl, header, units)
    !=======================================================================
    !        double precision
    !=======================================================================
    integer(I4b), parameter :: KCL = DP
    !
    CHARACTER(LEN=*),                          INTENT(IN) :: filename
    INTEGER(I4B),                              INTENT(IN) :: lmax, ncl
    REAL(KCL),         DIMENSION(0:lmax, 1:ncl), INTENT(OUT) :: clin
    CHARACTER(LEN=*), DIMENSION(1:),   INTENT(OUT) :: header
    CHARACTER(LEN=*), dimension(1:), optional,  INTENT(OUT) :: units

    print*,"-------------------------------------------------------------"
    print*,"WARNING : the routine read_asctab is now obsolete"
    print*,"  Use "
    print*," call fits2cl(filename, clin, lmax, ncl, header, units)"
    print*,"  instead (same module)"
    print*,"-------------------------------------------------------------"

    call fits2cl(filename, clin, lmax, ncl, header, units)
    return
  end subroutine read_asctab_d

  !=======================================================================
  subroutine output_map_s(map, header, outfile, extno)
    !=======================================================================
    ! this is only a wrapper to write_bintab for SP maps
    ! this routine needs an explicit interface
    !=======================================================================
    use long_intrinsic, only: long_size

    REAL(SP),         INTENT(IN), DIMENSION(0:,1:), target :: map
    CHARACTER(LEN=*), INTENT(IN), DIMENSION(1:) :: header
    CHARACTER(LEN=*), INTENT(IN)                :: outfile
    INTEGER(I4B)    , INTENT(IN)     , optional :: extno

    INTEGER(I8B) :: npix
    INTEGER(I4B) :: nlheader, nmap, extno_i
    real(SP), pointer, dimension(:,:) :: pmap
    !-----------------------------------------------------------------------
    npix = long_size(map,1)
    nmap = long_size(map,2)
    nlheader = SIZE(header)

    extno_i = 0
    if (present(extno)) extno_i = extno
    pmap => map(0:npix-1,1:nmap)
    call write_bintab(pmap, npix, nmap, &
         &            header(1:nlheader), nlheader, outfile, extno=extno_i)
    return
  end subroutine output_map_s
  !=======================================================================
  subroutine output_map_d(map, header, outfile, extno)
    !=======================================================================
    ! this is only a wrapper to write_bintab for DP maps
    ! this routine needs an explicit interface
    !=======================================================================
    use long_intrinsic, only: long_size

    REAL(DP),         INTENT(IN), DIMENSION(0:,1:), target :: map
    CHARACTER(LEN=*), INTENT(IN), DIMENSION(1:) :: header
    CHARACTER(LEN=*), INTENT(IN)                :: outfile
    INTEGER(I4B)    , INTENT(IN)     , optional :: extno

    INTEGER(I8B) :: npix
    INTEGER(I4B) :: nlheader, nmap, extno_i
    real(DP), pointer, dimension(:,:) :: pmap
    !-----------------------------------------------------------------------
    npix = long_size(map,1)
    nmap = long_size(map,2)
    nlheader = SIZE(header)

    extno_i = 0
    if (present(extno)) extno_i = extno
    pmap => map(0:npix-1,1:nmap)
    call write_bintab(pmap, npix, nmap, &
         &            header(1:nlheader), nlheader, outfile, extno=extno_i)
    return
  end subroutine output_map_d



  !=======================================================================
  subroutine write_dbintab(plm, nplm, nhar, header, nlheader, filename, nsmax, nlmax)
    !=======================================================================
    ! now obsolete routine
    !=======================================================================
    INTEGER(I4B), INTENT(IN) :: nplm, nhar, nlheader, nsmax, nlmax
    REAL(DP),          INTENT(IN), DIMENSION(0:nplm-1,1:nhar) :: plm
    CHARACTER(LEN=80), INTENT(IN), DIMENSION(1:nlheader) :: header
    CHARACTER(LEN=*),  INTENT(IN)               :: filename
    !=======================================================================
    print*,'WRITE_DBINTAB is obsolete.'
    print*,'   '
    print*,'To write a Healpix map into a FITS file'
    print*,'use WRITE_BINTAB or OUTPUT_MAP, which both accept '
    print*,'Single and Double Precision arguments'
    print*,'   '
    print*,'To write Plm polynoms into a FITS file,'
    print*,'use WRITE_PLM  (same syntax)'
    call write_plm(plm, nplm, nhar, header, nlheader, filename, nsmax, nlmax)

    return
  end subroutine write_dbintab
  !=======================================================================
  ! WRITE_PLM
  !     Create a FITS file containing a binary table extension with 
  !     the temperature map in the first column
  !     written by EH from writeimage and writebintable 
  !     (fitsio cookbook package)
  !
  !     slightly modified to deal with vector column 
  !     in binary table       EH/IAP/Jan-98
  !
  !     simplified the calling sequence, the header sould be filled in
  !     before calling the routine
  !
  !     Changed to write double precision plms for const.real. FKH/Apr-99
  !
  !=======================================================================
#ifndef NO64BITS
  subroutine write_plm4(plm, nplm, nhar, header, nlheader, filename, nsmax, nlmax)
    !=======================================================================
    INTEGER(I4B), INTENT(IN) :: nplm
    INTEGER(I4B), INTENT(IN) :: nhar, nlheader, nsmax, nlmax
    REAL(DP),          INTENT(IN), DIMENSION(0:nplm-1,1:nhar) :: plm
    CHARACTER(LEN=80), INTENT(IN), DIMENSION(1:nlheader) :: header
    CHARACTER(LEN=*),  INTENT(IN)               :: filename
    integer(i8b) :: nplm8
    nplm8 = nplm
    call write_plm8(plm, nplm8, nhar, header, nlheader, filename, nsmax, nlmax)
    return
  end subroutine write_plm4
    !=======================================================================
#endif
  subroutine write_plm8(plm, nplm, nhar, header, nlheader, filename, nsmax, nlmax)
    !=======================================================================

    INTEGER(I8B), INTENT(IN) :: nplm
    INTEGER(I4B), INTENT(IN) :: nhar, nlheader, nsmax, nlmax
    REAL(DP),          INTENT(IN), DIMENSION(0:nplm-1,1:nhar) :: plm
    CHARACTER(LEN=80), INTENT(IN), DIMENSION(1:nlheader) :: header
    CHARACTER(LEN=*),  INTENT(IN)               :: filename
    !
    INTEGER(I4B) ::  status,unit,blocksize,bitpix,naxis,naxes(1)
    INTEGER(I4B) ::  i
    LOGICAL(LGT) ::  simple,extend
    CHARACTER(LEN=80) :: comment

    INTEGER(I4B), PARAMETER :: MAXDIM = MAXDIM_TOP !number of columns in the extension
    INTEGER(I4B) :: nrows, tfields, varidat, nmmax
    INTEGER(I8B) :: nlm, nf, i0, i1
    INTEGER(I4B) :: frow,  felem, colnum
    CHARACTER(LEN=20) :: ttype(MAXDIM), tform(MAXDIM), tunit(MAXDIM), extname
    CHARACTER(LEN=10) ::  card, num
    CHARACTER(LEN=2)  :: stn
    INTEGER(I4B)      :: itn, irow, np
    real(DP)          :: dm
    !-----------------------------------------------------------------------

    status=0

    unit = 148

    !     create the new empty FITS file
    blocksize=1
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

    nlm = nplm / (2*nsmax) ! number of Plm per Healpix ring
    dm  = (2*nlmax + 3.0_dp)**2 - 8.0_dp*nlm
    nmmax = nlmax - (nint(sqrt(dm)) - 1)/2
    call assert(nplm == nlm*2*nsmax, 'inconsistent array size in write_plm')
    call assert(nplm == (nmmax+1_i8b)*(2*nlmax-nmmax+2_i8b)*nsmax, &
         &     'inconsistent array size (nmmax) in write_plm')

    !     creates an extension
    call ftcrhd(unit, status)

    !     writes required keywords
    nrows    = 2*nsmax ! naxis1
    tfields  = nhar
    nf = nplm / nrows
    if (nf * nrows /= nplm) then
       print*,' problems in write_plm',nf*nrows,nplm
       stop
    endif
    write (num,'(I10)') nf
!     tform(1:nhar) = trim(adjustl(num))//'D'     ! does not work with Ifort, EH, 2006-04-04
    num = trim(adjustl(num))//'D'
    tform(1:nhar) = num
    ttype(1:nhar) = 'harmonics'   ! will be updated
    tunit(1:nhar) = ''      ! optional, will not appear
    extname  = ''      ! optional, will not appear
    varidat  = 0
    call ftphbn(unit, nrows, tfields, ttype, tform, tunit, &
         &     extname, varidat, status)

    !     write the header literally, putting TFORM1 at the desired place
    comment = 'data format of field: 8-byte REAL'
    do i=1,nlheader
       card = header(i)
       if (card(1:5) == 'TTYPE') then ! if TTYPE1 is explicitely given
          stn = card(6:6)
          read(stn,'(i1)') itn
          ! discard at their original location:
          call ftdkey(unit,'TTYPE'//stn,status)  ! old TTYPEi and
          status = 0
          call ftdkey(unit,'TFORM'//stn,status)  !     TFORMi
          status = 0
          if (itn <= tfields) then ! only put relevant information 2008-08-27
             call putrec(unit,header(i), status)           ! write new TTYPE1
             status = 0
             call ftpkys(unit,'TFORM'//stn,tform(itn),comment,status) ! and write new TFORM1 right after
          endif
       elseif (header(i)/=' ') then
          call putrec(unit,header(i), status)
       endif
       status = 0
    enddo

    ! make sure keywords critical for future use are present and correct
    call ftukyj(unit,'NSIDE',   nsmax, 'HEALPIX resolution parameter', status)
    call ftukyj(unit,'MAX-LPOL',nlmax, 'Maximum multipole order l of P_lm', status)
    call ftukyj(unit,'MAX-MPOL',nmmax, 'Maximum degree m of P_lm', status)
    call ftukyl(unit,'POLAR', (nhar>1),'Polarisation included (T/F)', status)

    !     write the extension one column by one column
    do irow=1, nrows
       frow   = irow  ! starting position (row)
       felem  = 1  ! starting position (element)
       np = nf
       i0 = (irow-1)*nf
       i1 = i0+nf-1
       do colnum = 1, nhar
          !call ftpcld(unit, colnum, frow, felem, nplm, plm(0,colnum), status)
          call ftpcld(unit, colnum, frow, felem, np, plm(i0:i1,colnum), status)
       enddo
    enddo

    !     ----------------------
    !     close and exit
    !     ----------------------

    !     close the file and free the unit number
    call ftclos(unit, status)

    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)

    return
  end subroutine write_plm8

  !=======================================================================
  !       ALMS2FITS
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
  !=======================================================================
  subroutine alms2fits_s(filename, nalms, alms, ncl, header, nlheader, next)
    !=======================================================================
    character(len=*),  intent(in)               :: filename
    integer(I4B),      intent(in)               :: nalms, nlheader, next, ncl
    real(SP),          intent(in), dimension(1:nalms,1:ncl+1,1:next), target :: alms
    character(len=80), intent(in), dimension(1:nlheader,1:next) :: header
    integer(I4B) :: i
    real(SP),          dimension(:,:), pointer :: palms

    do i=1,next
       palms => alms(1:nalms,1:(ncl+1),i) ! use pointer to reduce memory usage
       call write_alms_s(filename, nalms, palms, ncl, &
            &            header(1:nlheader,i), nlheader, i)
    enddo

  end subroutine alms2fits_s
  !=======================================================================
  subroutine alms2fits_d(filename, nalms, alms, ncl, header, nlheader, next)
    !=======================================================================
    character(len=*),  intent(in)               :: filename
    integer(I4B),      intent(in)               :: nalms, nlheader, next, ncl
    real(DP),          intent(in), dimension(1:nalms,1:ncl+1,1:next), target :: alms
    character(len=80), intent(in), dimension(1:nlheader,1:next) :: header
    integer(I4B) :: i
    real(DP),          dimension(:,:), pointer :: palms

    do i=1,next
       palms => alms(1:nalms,1:(ncl+1),i) ! use pointer to reduce memory usage
       call write_alms_d(filename, nalms, palms, ncl, &
            &            header(1:nlheader,i), nlheader, i)
    enddo

  end subroutine alms2fits_d

  !=======================================================================
  subroutine fits2alms_s(filename, nalms, alms, ncl, header, nlheader, next) 
    !=======================================================================
    character(len=*),  intent(in)               :: filename
    integer(I4B),      intent(in)               :: nalms, nlheader, next, ncl
    real(SP),          intent(inout), dimension(1:nalms,1:ncl+1,1:next), target :: alms
    character(len=80), intent(out), dimension(1:nlheader,1:next) :: header
    integer(I4B) :: i
    real(SP),          dimension(:,:), pointer :: palms

    do i=1,next
       palms => alms(1:nalms,1:(ncl+1),i) ! use pointer to reduce memory usage
       call read_alms_s(filename, nalms, palms, ncl, &
            &           header(1:nlheader,i), nlheader, i)
    enddo

  end subroutine fits2alms_s
  !=======================================================================
  subroutine fits2alms_d(filename, nalms, alms, ncl, header, nlheader, next) 
    !=======================================================================
    character(len=*),  intent(in)               :: filename
    integer(I4B),      intent(in)               :: nalms, nlheader, next, ncl
    real(DP),          intent(inout), dimension(1:nalms,1:ncl+1,1:next), target :: alms
    character(len=80), intent(out), dimension(1:nlheader,1:next) :: header
    integer(I4B) :: i
    real(DP),          dimension(:,:), pointer :: palms

    do i=1,next
       palms => alms(1:nalms,1:(ncl+1),i) ! use pointer to reduce memory usage
       call read_alms_d(filename, nalms, palms, ncl, &
            &           header(1:nlheader,i), nlheader, i)
    enddo

  end subroutine fits2alms_d

  SUBROUTINE input_tod_s(filename, tod, npixtot, ntods, &
       &                 header, firstpix, fmissval, extno)
  ! *************************************************************************

  ! *************************************************************************
  ! v1.0 : Olivier Dore et Eric Hivon, June-July 2002
  ! v1.1     2002-07-08 : bugs correction by E.H. 
  ! v1.2     2002-08-22 : pixel number starts at 0 for these routine,
  !                         but starts at 1 for cftisio routines ...
  ! *************************************************************************
    ! ===================================================================================
    !     reads fits file
    !     filename = fits file (input)
    !     tod      = data rad from the file (ouput) = real*4 array of size (npixtot,ntods)
    !     npixtot  = length of the tod (input)
    !     ntods    = number of tods
    !     header   = OPTIONAL, FITS header
    !     fmissval = OPTIONAL argument (input) with the value to be given to missing
    !             pixels, its default value is 0
    !     firstpix = OPTIONAL first pixel to be read (starting at 0)
    !     extno    = OPTIONAL, extension number (starting at 0)
    !                                               O.D. & E.H. 06/02 @ IAP
    !      2002-07-08 : bugs correction by E.H. 
    !       (consistent use of fmiss_effct)
    ! =======================================================================
    
    IMPLICIT NONE

    INTEGER(I8B),     INTENT(IN)           :: npixtot
    INTEGER(I8B),     INTENT(IN), OPTIONAL :: firstpix
    INTEGER(I4B),     INTENT(IN)           :: ntods
    REAL(SP),         INTENT(OUT)          ,dimension(0:,1:) :: tod
    CHARACTER(LEN=*), INTENT(IN)           :: filename
    REAL(SP),         INTENT(IN), OPTIONAL :: fmissval
    character(len=*), intent(out),OPTIONAL, dimension(1:)  :: header
    INTEGER(I4B),     INTENT(IN), OPTIONAL :: extno
    
    INTEGER(I4B) :: i,itod
    REAL(SP)     :: fmissing, fmiss_effct
    INTEGER(I8B) :: imissing
    integer(i8b) :: fp = 0_i8b

    LOGICAL(LGT) :: anynull
    
    !-----------------------------------------------------------------------
    fmiss_effct = 0.
    IF (PRESENT(fmissval)) fmiss_effct = fmissval

    if (present(firstpix)) fp = firstpix

    CALL read_bintod_s(filename, tod, npixtot, ntods, fp, fmissing, anynull, &
         &             header, extno)

    call map_bad_pixels(tod, fmissing, fmiss_effct, imissing, verbose=.true.)
!     DO itod = 1, ntods
!        anynull = .TRUE.
!        IF (anynull) THEN
!           imissing = 0
!           DO i=0,npixtot-1
!              IF ( ABS(tod(i,itod)/fmissing -1.) .LT. 1.e-5 ) THEN
!                 tod(i,itod) = fmiss_effct
!                 imissing = imissing + 1
!              ENDIF
!           ENDDO
!           IF (imissing .GT. 0) THEN
!              WRITE(*,'(a,1pe11.4)') 'blank value : ' ,fmissing
!              WRITE(*,'(i7,a,f7.3,a,1pe11.4)') &
!                   &           imissing,' missing pixels (', &
!                   &           (100.*imissing)/npixtot,' %),'// &
!                   &           ' have been set to : ',fmiss_effct
!           ENDIF
!        ENDIF
!     ENDDO
    RETURN

  END SUBROUTINE input_tod_s
  !=======================================================================

  !**************************************************************************************
  SUBROUTINE input_tod_d(filename, tod, npixtot, ntods, &
       &                 header, firstpix, fmissval, extno)
  !**************************************************************************************

    !=======================================================================
    !     reads fits file
    !     filename = fits file (input)
    !     tod      = data rad from the file (ouput) = real*8 array of size (npixtot,ntods)
    !     npixtot  = length of the tod (input)
    !     ntods     = number of tods
    !     header   = OPTIONAL, FITS header
    !     fmissval = OPTIONAL argument (input) with the value to be given to missing
    !             pixels, its default value is 0
    !     firstpix = OPTIONAL first pixel to be read (starting at 0)
    !     extno    = OPTIONAL, extension number (starting at 0)
    !                                               O.D. & E.H. 06/02 @ IAP
    !      2002-07-08 : bugs correction by E.H. 
    !       (consistent use of fmiss_effct)
    !=======================================================================
    
    IMPLICIT NONE

    INTEGER(I8B),     INTENT(IN)           :: npixtot
    INTEGER(I8B),     INTENT(IN), OPTIONAL :: firstpix
    INTEGER(I4B),     INTENT(IN)           :: ntods
    REAL(DP),         INTENT(OUT)          ,dimension(0:,1:) :: tod
    CHARACTER(LEN=*), INTENT(IN)           :: filename
    REAL(DP),         INTENT(IN), OPTIONAL :: fmissval
    character(len=*), intent(out),OPTIONAL, dimension(1:)  :: header
    INTEGER(I4B),     INTENT(IN), OPTIONAL :: extno
    
    INTEGER(I4B) :: i,itod
    REAL(DP)     :: fmissing, fmiss_effct
    INTEGER(I8B) :: imissing
    integer(i8b) :: fp = 0_i8b

    LOGICAL(LGT) :: anynull
    
    !-----------------------------------------------------------------------
    
    fmiss_effct = 0.
    IF (PRESENT(fmissval)) fmiss_effct = fmissval

    if (present(firstpix)) fp = firstpix
    
    CALL read_bintod_d(filename, tod, npixtot, ntods, fp, fmissing, anynull, &
         &             header, extno)

    call map_bad_pixels(tod, fmissing, fmiss_effct, imissing, verbose=.true.)
!     DO itod = 1, ntods
!        anynull = .TRUE.
!        IF (anynull) THEN
!           imissing = 0
!           DO i=0,npixtot-1
!              IF ( ABS(tod(i,itod)/fmissing -1.) .LT. 1.e-5 ) THEN
!                 tod(i,itod) = fmiss_effct
!                 imissing = imissing + 1
!              ENDIF
!           ENDDO
!           IF (imissing .GT. 0) THEN
!              WRITE(*,'(a,1pe11.4)') 'blank value : ' ,fmissing
!              WRITE(*,'(i7,a,f7.3,a,1pe11.4)') &
!                   &           imissing,' missing pixels (', &
!                   &           (100.*imissing)/npixtot,' %),'// &
!                   &           ' have been set to : ',fmiss_effct
!           ENDIF
!        ENDIF
!     ENDDO
    RETURN

  END SUBROUTINE input_tod_d
  !=======================================================================
  !=======================================================================
  subroutine read_dbintab(filename,map,npixtot,nmaps,nullval,anynull, units)
    !=======================================================================
    !     Read a FITS file
    !
    !     slightly modified to deal with vector column 
    !     in binary table       EH/IAP/Jan-98
    !
    !     Reads a double-precision array with precomputed plms, used by syn/anafast
    !                FKH/Apr-99
    !=======================================================================
    CHARACTER(LEN=*),                           INTENT(IN) :: filename
    INTEGER(I4B),                               INTENT(IN) :: npixtot, nmaps
    !       REAL(DP), DIMENSION(0:npixtot-1,1:nmaps),   INTENT(OUT) :: map
    REAL(DP), DIMENSION(0:,1:),                 INTENT(OUT) :: map
    REAL(DP),                                   INTENT(OUT) :: nullval
    LOGICAL(LGT),                               INTENT(OUT) ::  anynull
    !       CHARACTER(LEN=*), dimension(1:), optional,           INTENT(OUT) :: header
    CHARACTER(LEN=*), dimension(1:), optional,  INTENT(OUT) :: units

    INTEGER(I4B) :: status,unit,readwrite,blocksize,naxes(2),nfound, naxis
    INTEGER(I4B) :: group,firstpix,npix,i
    REAL(SP) :: blank, testval
    REAL(DP) ::  bscale,bzero
    CHARACTER(LEN=80) :: comment
    LOGICAL(LGT) :: extend
    INTEGER(I4B) :: nmove, hdutype
    INTEGER(I4B) :: column, frow, imap
    INTEGER(I4B) :: datacode, repeat, width

    INTEGER(I4B), PARAMETER :: MAXDIM = MAXDIM_TOP !number of columns in the extension
    INTEGER(I4B) :: nrows, tfields, varidat
    CHARACTER(LEN=20), dimension(1:MAXDIM) :: ttype, tform, tunit
    CHARACTER(LEN=20)                      :: extname
    !-----------------------------------------------------------------------
    status=0

    unit = 147
    naxes(1) = 1
    naxes(2) = 1
    nfound = -1
    anynull = .false.
    bscale = 1.0d0
    bzero = 0.0d0
    blank = -2.e25
    nullval = bscale*blank + bzero

    readwrite=0
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

    if (naxis > 0) then ! there is an image
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

       npix=naxes(1)
       if (npix /= npixtot) then
          print *,'found ',npix,' plms'
          print *,'expected ',npixtot
          call fatal_error
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
       call ftgkye(unit,'BLANK', blank, comment,status)
       if (status == 202) then ! BLANK not found 
          ! (according to fitsio BLANK is integer)
          blank = -2.e25
          status = 0
       endif
       nullval = bscale*blank + bzero

       !        -----------------------------------------

       group=1
       firstpix = 1
       call ftgpvd(unit,group,firstpix,npix,nullval,map,anynull,status)
       ! if there are any NaN pixels, (real data)
       ! or BLANK pixels (integer data) they will take nullval value
       ! and anynull will switch to .true.
       ! otherwise, switch it by hand if necessary
       testval = 1.e-6 * ABS(nullval)
       do i=0, npix-1
          if (ABS(map(i,1)-nullval) < testval) then
             anynull = .true.
             goto 111
          endif
       enddo
111    continue

    else if (extend) then ! there is an extension
       nmove = +1
       call ftmrhd(unit, nmove, hdutype, status)
       !cc         write(*,*) hdutype

       call assert (hdutype==2, 'this is not a binary table')

       !        reads all the keywords
       call ftghbn(unit, MAXDIM, &
            &        nrows, tfields, ttype, tform, tunit, extname, varidat, &
            &        status)

       if (tfields < nmaps) then
          print *,'found ',tfields,' maps in the file'
          print *,'expected ',nmaps
          call fatal_error
       endif

       !        finds the bad data value
       call ftgkyd(unit,'BAD_DATA',nullval,comment,status)
       if (status == 202) then ! bad_data not found
          nullval = d_bad_value
          status = 0
       endif

       !          if (nlheader > 0) then
       !             call get_clean_header(unit, header, filename, status)
       !          endif

       if (present(units)) then
          units = 'unknown'
          do imap = 1, nmaps
             units(imap) = adjustl(tunit(imap))
          enddo
       endif

       do imap = 1, nmaps
          !parse TFORM keyword to find out the length of the column vector
          call ftbnfm(tform(imap), datacode, repeat, width, status)

          !reads the columns
          column = imap
          frow = 1
          firstpix = 1
          npix = nrows * repeat
          if (npix /= npixtot) then
             print *,'found ',npix,' plms'
             print *,'expected ',npixtot
             call fatal_error("read_dbintab "//trim(filename))
          endif
          call ftgcvd(unit, column, frow, firstpix, npix, nullval, &
               &        map(0:npix-1,imap), anynull, status)
       enddo

    else ! no image no extension, you are dead, man
       call fatal_error(' No image, no extension')
    endif

    !     close the file
    call ftclos(unit, status)

    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)

    return
  end subroutine read_dbintab


  !=======================================================================
  subroutine printerror(status)
    !=======================================================================
    !     Print out the FITSIO error messages to the user
    !=======================================================================
    INTEGER(I4B), INTENT(IN) :: status
    CHARACTER ::  errtext*30,errmessage*80
    !-----------------------------------------------------------------------
    !     check if status is OK (no error); if so, simply return
    if (status .le. 0)return

    !     get the text string which describes the error
    call ftgerr(status,errtext)
    print *,'FITSIO Error Status =',status,': ',errtext

    !     read and print out all the error messages on the FITSIO stack
    call ftgmsg(errmessage)
    do while (errmessage /= ' ')
       print *,errmessage
       call ftgmsg(errmessage)
    end do

    return
  end subroutine printerror


  !=======================================================================
  subroutine read_par(filename,nside,lmax,tfields,mmax)
    !=======================================================================
    !        Read nside, lmax, tfields and mmax from a FITS file          
    !    parameters not found take a value of -1
    !
    !         Frode K. Hansen, April 1999
    !         EH, Dec 2004
    !
    !=======================================================================
    CHARACTER(LEN=*), INTENT(IN)            :: filename

    INTEGER(I4B), INTENT(OUT)   :: nside, lmax, tfields
    integer(i4b), intent(out), optional :: mmax

    INTEGER(I4B) :: status,unit,readwrite,blocksize,naxis
    CHARACTER(LEN=80) :: comment, ttype1
    LOGICAL(LGT) ::  extend, anyf
    INTEGER(I4B)::  nmove, hdutype, nrows
    integer(I4B), dimension(1) :: idmax ! because ftgcvj expects an array (2020-08-24)

    !-----------------------------------------------------------------------
    status=0
    unit = 120
    comment=''
    ttype1=''

    readwrite=0
    call ftopen(unit,filename,readwrite,blocksize,status)
    if (status > 0) call printerror(status)
    !     -----------------------------------------
    call ftgkyj(unit,'NAXIS', naxis, comment, status)

    !     determines the presence of an extension
    call ftgkyl(unit,'EXTEND', extend, comment, status)
    call assert (status<=0, 'No Extension in FITS file!')

    nmove = +1
    call ftmrhd(unit, nmove, hdutype, status)

    call assert (hdutype==2, 'This is not a FITS binary-table')

    call ftgkyj(unit,'NSIDE',nside,comment,status)
    if (status == 202) then
       print*,'WARNING: NSIDE keyword not found!'
       nside = -1
       status = 0
    endif

    call ftgkyj(unit,'TFIELDS',tfields,comment,status)
    if (status == 202) then
       print*,'WARNING: TFIELDS keyword not found!'
       tfields = -1
       status = 0
    endif
    
    call ftgkyj(unit,'MAX-LPOL',lmax,comment,status)
    if (status == 202) then
       status = 0
       ! if not found, determines if file contains indexed list of alm
       ! if so, find out largest l
       if (tfields >= 3 .and. hdutype==2) then ! 3 column binary table
          call ftgkys(unit,'TTYPE1',ttype1,comment,status)
          ttype1 = trim(strupcase(adjustl(ttype1)))
          if (trim(ttype1(1:5)) == 'INDEX') then
             call ftgkyj(unit, 'NAXIS2', nrows, comment, status) ! find number of rows
             call ftgcvj(unit, 1_i4b, nrows, 1_i4b, 1_i4b, 0_i4b, idmax(1), anyf, status) ! read element on last row of first column
             if (status == 0) then
                lmax = int(sqrt(   real(idmax(1)-1, kind = DP)  ) )
                if (lmax > 0) goto 1000
             endif
          endif
       endif
       print*,'WARNING: MAX-LPOL keyword not found!'
       lmax = -1
       status = 0
    endif
1000 continue

    if (present(mmax)) then
       call ftgkyj(unit,'MAX-MPOL',mmax,comment,status)
       if (status == 202) then
          print*,'WARNING: MAX-MPOL keyword not found!'
          mmax = -1
          status = 0
       endif
    endif
    call ftclos(unit, status)

  end subroutine read_par

  !=======================================================================
  function getnumext_fits(filename)
    !=======================================================================
    !  result = getnumext_fits(filename)
    !    returns the number of extensions present in FITS file 'filename'
    !
    ! EH, Nov 2004
    ! April 2007: close file on exit
    !=======================================================================
    character(LEN=*), intent(IN)             :: filename
    integer(i4b)                             :: getnumext_fits
    !
    integer(i4b) :: status, unit, readwrite, blocksize, nhdu
    !-----------------------------------------------------------------------
    status         = 0
    unit           = 149
    getnumext_fits = 0
    readwrite      = 0 ! Read only
    call ftopen(unit, filename, readwrite, blocksize, status)
    if (status > 0) then
       call printerror(status)
       call ftclos(unit, status)
       return
    endif

    call ftthdu(unit, nhdu, status)
    getnumext_fits = nhdu - 1

    call ftclos(unit, status)
    return
  end function getnumext_fits
  !=======================================================================
  function getsize_fits(filename, nmaps, ordering, obs_npix, &
       &    nside, mlpol, type, polarisation, fwhm_arcmin, beam_leg, & 
       &    coordsys, polcconv, extno)
    !=======================================================================
    !  result = getsize_fits(filename, nmaps, ordering, nside, mlpol, type)
    !     Get the number of pixels stored in a map FITS file.
    !     Each pixel is counted only once 
    !     (even if several information is stored on each of them, see nmaps).
    !     Depending on the data storage format, this may be :
    !       - equal or smaller to the number Npix of Healpix pixels available 
    !          over the sky for the given resolution (Npix = 12*nside*nside)
    !       - equal or larger to the number of non blank pixels (obs_npix)
    !     
    !     filename = (IN) name of the FITS file containing the map
    !     nmaps = (OPTIONAL, OUT) number of maps in the file
    !     ordering = (OPTIONAL, OUT) Healpix ordering scheme used
    !                  0 = unknown
    !                  1 = RING
    !                  2 = NESTED
    !     obs_npix =  (OPTIONAL, OUT) number of non blanck pixels.    
    !                It is set to -1 if it can not be determined from header
    !                information alone
    !     nside = (OPTIONAL, OUT) Healpix parameter Nside
    !                 returns a negative value if not found
    !     mlpol = (OPTIONAL, OUT) maximum multipole used to generate the map (for simulated map)
    !                 returns a negative value if not found
    !     type = (OPTIONAL, OUT) Healpix/FITS file type
    !                  <0 : file not found, or not valid
    !                  0  : image only fits file, deprecated Healpix format
    !                        (result = 12 * nside * nside)
    !                  1  : ascii table, generally used for C(l) storage
    !                  2  : binary table : with implicit pixel indexing (full sky)
    !                        (result = 12 * nside * nside)
    !                  3  : binary table : with explicit pixel indexing (generally cut sky)
    !                        (result <= 12 * nside * nside)
    !                999  : unable to determine the type
    !     polarisation = (OPTIONAL, OUT) presence of polarisation data in the file
    !                  <0 : can not find out
    !                   0 : no polarisation
    !                   1 : contains polarisation (Q,U or G,C)
    !     fwhm_arcmin     = (OPTIONAL, DP, OUT) returns the beam FWHM read from FITS header, 
    !                        translated from Deg (hopefully) to arcmin
    !                     returns a negative value if not found
    !
    !     beam_leg     = (OPTIONAL, CHR, OUT) filename of beam or filtering window function applied to data
    !                     returns a empty string if not found
    !     coordsys     = (OPTIONAL, CHR, OUT) string describing coordinate system,  
    !                     G = Galactic, E = ecliptic, C = celestial = equatorial.
    !                     empty if not found.
    !
    !     polcconv     = (OPTIONAL, I4B, OUT) coordinate convention for polarisation
    !                   0: not found
    !                   1: COSMO (default for Healpix)
    !                   2: IAU
    !                   3: neither COSMO nor IAU
    !
    !     extno = (OPTIONAL, IN) specify FITS extension to look at (0 based)
    !                  default = 0 (first extension)
    !
    !     Benjamin D. Wandelt January 1998
    !     includes Eric Hivon's modification of FITS columns
    !     and draws heavily on the read_bintab routine.
    !
    !     addition of optional Ordering by E.H. (July 98)
    !     addition of optional Nside and Mlpol by ?? (??)
    !     addition of optional type by E.H. (Sept 00)
    !     improved for ASCII table E.H (Feb 03)
    !     addition of extno, E.H (Nov 04)
    !     addition of fwhm_arcmin and beam_leg, EH (Jan 05).
    !     addition of polcconv, EH (June 05).
    !     accept files with NAXIS2 > 2^31 (2015-04) (requires recent cfitsio >= 3.20)
    !     accept files with polarisation columns named Q_S*, Q-S*, ... U_S*, U-S*
    !=======================================================================
    use pix_tools, only: nside2npix
    character(LEN=*), intent(IN)             :: filename
    integer(I4B),     intent(out),  optional :: nmaps
    integer(I4B),     intent(out),  optional :: ordering
    integer(I4B),     intent(out),  optional :: nside
    integer(I4B),     intent(out),  optional :: mlpol
    integer(I4B),     intent(out),  optional :: obs_npix
    integer(I4B),     intent(out),  optional :: type
    integer(I4B),     intent(out),  optional :: polarisation
    real(DP),         intent(out),  optional :: fwhm_arcmin
    character(LEN=*), intent(out),  optional :: beam_leg
    character(LEN=*), intent(out),  optional :: coordsys
    integer(I4B),     intent(out),  optional :: polcconv
    integer(I4B),     intent(in),   optional :: extno

    INTEGER(I4B)           :: nmaps_in, ordering_in, nside_in
    INTEGER(I4B)           :: mlpol_in, obs_npix_in, ftype_in
    INTEGER(I4B)           :: extno_in, polcconv_in
    real(DP)               :: fwhm_arcmin_in
    character(len=FILENAMELEN)        :: beam_leg_in
    character(len=20)        :: coordsys_in
    INTEGER(I4B)           :: grain
    CHARACTER(len=20)      :: indxschm
    CHARACTER(LEN=20)      :: order_val, object_val, ttype_val, polcconv_val
!     INTEGER(I4B)           :: getsize_fits
    INTEGER(I8B)           :: getsize_fits
    LOGICAL(LGT)           ::  polar_in
    character(len=3),  dimension(1:16,1:2)  :: defpol
    logical(kind=LGT), dimension(1:2)       :: pf
    integer(kind=I4B)                       :: ndp, j, k

    INTEGER(I4B)      :: status,unit,readwrite,blocksize,naxes(2),nfound, naxis
    INTEGER(I4B)      :: i
!     INTEGER(I4B) :: nsize
    INTEGER(I8B)      :: nsize
    CHARACTER(LEN=80) :: comment
    LOGICAL(LGT)      ::  extend
    INTEGER(I4B)      ::  nmove, hdutype, hdunum
    INTEGER(I4B)      :: datacode, repeat1, repeat2, width

    INTEGER(I4B),           PARAMETER :: MAXDIM = MAXDIM_TOP !number of columns in the extension
    INTEGER(I4B)                      :: tfields, varidat, rowlen
    CHARACTER(LEN=20), dimension(1:MAXDIM) :: ttype, tform, tunit
    INTEGER(I4B),      dimension(1:MAXDIM) :: tbcol
    CHARACTER(LEN=20)                      :: extname
    integer(I4B) :: nrows32
    integer(I8B) :: nrows64
    !-----------------------------------------------------------------------
    status=0
    order_val = ''
    nmaps_in = 0
    ordering_in = 0
    mlpol_in = -1
    nside_in = -1
    ftype_in = 999 !   
    obs_npix_in = -1
    unit = 119
    naxes(1) = 1
    naxes(2) = 1
    nfound = -1
    extno_in = 0
    polcconv_in = 0
    comment=''
    ttype=''
    tform=''
    tunit=''
    extname=''
    if (present(extno)) extno_in = extno
    naxis = 0
    extend = .false.

    readwrite=0
    !call ftopen(unit,filename,readwrite,blocksize,status)
    call ftnopn(unit,filename,readwrite,status)
    if (status > 0) then
       ftype_in = -1
       getsize_fits = -1
       call printerror(status)
       call ftclos(unit, status)
       return
    endif
    !     -----------------------------------------
    call ftghdn(unit, hdunum)
    if (hdunum == 1) then  ! in primary HDU: move to next HDU
       !     determines the presence of image
       call ftgkyj(unit,'NAXIS', naxis, comment, status)

       !     determines the presence of an extension
       call ftgkyl(unit,'EXTEND', extend, comment, status)
       if (status > 0) then
          ftype_in = 0
          status = 0 ! no extension : 
       !     to be compatible with first version of the code
       endif

    else ! already in non primary HDU
       extend = .true.
       call ftghdt(unit, hdutype, status)
    endif

    if (naxis > 0) then 
       !---------------------------------
       ! there is an image
       !---------------------------------
       !        determine the size of the image (look naxis1 and naxis2)
       call ftgknj(unit,'NAXIS',1_i4b,2_i4b,naxes,nfound,status)

       !        check that it found only NAXIS1
       if (nfound == 2 .and. naxes(2) > 1) then
          print *,'multi-dimensional image'
          print *,'expected 1-D data.'
          call ftclos(unit, status)
          call fatal_error
       end if

       if (nfound < 1) then
          call printerror(status)
          print *,'can not find NAXIS1.'
          call ftclos(unit, status)
          call fatal_error
       endif

       nsize=naxes(1)

       call ftgkys(unit,'ORDERING',order_val,comment,status)
       if (status == 202) then ! Ordering not found
          ordering_in = 0
          order_val = ''
          status = 0
       endif

       call ftgkyj(unit,'NSIDE',nside_in,comment,status)
       if (status == 202) then ! Nside not found
          nside_in = -1
          status = 0
       endif

       if (present(mlpol)) then
          call ftgkyj(unit,'MAX-LPOL',mlpol_in,comment,status)
          if (status == 202) then ! max-lpol not found
             mlpol_in = -1
             status = 0
          endif
       endif

       if (present(polarisation)) then
          polarisation = 0
       endif
    else if (extend) then 
       !-----------------------------------------------------------------
       ! there is an extension
       !-----------------------------------------------------------------

       nmove =  extno_in
       if (hdunum == 1) nmove = extno_in + 1
       call ftmrhd(unit, nmove, hdutype, status)
       if (status > 0) then ! extension not found
          print*,'Extension #',extno_in,' not found in '//trim(filename)
          call printerror(status)
          call ftclos(unit, status)
          call fatal_error
       endif
       !c         write(*,*) hdutype

       !        reads all the keywords
       if (hdutype == 2) then ! binary table
          call ftghbn(unit, MAXDIM, &
            &         nrows32, tfields, ttype,        tform, tunit, &
            &         extname, varidat, status)
       else ! ASCII table (hdutype = 1)
          ftype_in = 1
          call ftghtb(unit, MAXDIM, &
            &         rowlen, &
            &         nrows32, tfields, ttype, tbcol, tform, tunit, &
            &         extname,          status)
       endif
       ! deal with NAXIS2 > 2^31, 2015-04
#ifdef NO64BITS
       call ftgnrw  (unit, nrows32, status)
       nrows64 = nrows32
#else
       call ftgnrwll(unit, nrows64, status)
#endif

       !        parse TFORM keyword to find out the length of the column vector
       repeat1 = 1
       repeat2 = 1
       call ftbnfm(tform(1), datacode, repeat1, width, status)
       if (tfields > 1) call ftbnfm(tform(2), datacode, repeat2, width, status)

       !nsize = int(nrows32, kind=i8b) * max(repeat1,repeat2) ! corrected Oct-03
       nsize = int(nrows64, kind=i8b) * max(repeat1,repeat2) ! 2015-04

       nmaps_in = tfields

       call ftgkys(unit,'ORDERING',order_val,comment,status)
       if (status == 202) then ! Ordering not found
          ordering_in = 0
          order_val = ''
          status = 0
       endif

       call ftgkyj(unit,'NSIDE',nside_in,comment,status)
       if (status == 202) then ! Nside not found
          nside_in = -1
          status = 0
       endif

       call ftgkyj(unit,'OBS_NPIX',obs_npix_in,comment,status)
       if (status == 202) then ! obs_npix not found
          obs_npix_in = -1
          status = 0
       endif

       if (present(mlpol)) then
          call ftgkyj(unit,'MAX-LPOL',mlpol_in,comment,status)
          if (status == 202) then ! max-lpol not found
             status = 0
             call ftgkyj(unit,'LMAX',mlpol_in,comment,status)
             if (status == 202) then
                mlpol_in = -1
                status = 0
             endif
          endif
       endif

       if (present(fwhm_arcmin)) then
          call ftgkyd(unit,'FWHM',fwhm_arcmin_in, comment, status)
          if (status == 202) then ! fwhm not found
             fwhm_arcmin_in = -1.
             status = 0
          else
             fwhm_arcmin_in = 60.0_dp * fwhm_arcmin_in
          endif
       endif

       if (present(beam_leg)) then
          call ftgkys(unit,'BEAM_LEG',beam_leg_in, comment, status)
          if (status == 202) then ! beam_leg not found
             beam_leg_in = ' '
             status = 0
          endif
       endif

       if (present(coordsys)) then
          call ftgkys(unit,'COORDSYS',coordsys_in, comment, status)
          if (status == 202) then ! coordsys not found
             coordsys_in = ' '
             status = 0
          endif
       endif

       ! determines pixel indexing (for binary tables)
       if (present(type) .and. ftype_in == 999) then
          ! most stringent test
          call ftgkys(unit,'INDXSCHM',indxschm,comment,status)
          if (status == 0) then ! found
             ftype_in = 3
             if (trim(indxschm) == 'IMPLICIT') ftype_in = 2
             goto 1000
          else
             status = 0
          endif
          ! 2nd most stringent test
          call ftgkyj(unit,'GRAIN',grain,comment,status)
          if (status == 0) then ! found
             ftype_in = 3
             if (grain == 0) ftype_in = 2
             goto 1000
          else
             status = 0
          endif
          ! 3rd most stringent test
          if (trim(ttype(1)) /= '') then
             if (trim(ttype(1)) == 'PIXEL' .or. trim(ttype(1)) == 'PIX') then
                ftype_in = 3
             else
                ftype_in = 2
             endif
             goto 1000
          endif
          ! lousy test
          call ftgkys(unit,'OBJECT',object_val,comment,status)
          if (status == 0) then
             if (trim(object_val) == 'PARTIAL') ftype_in = 3
             if (trim(object_val) == 'FULLSKY') ftype_in = 2
             if (ftype_in /= 999) goto 1000
          else
             status = 0
          endif
          ! very lousy test
          if (nside_in > 0) then
             ftype_in = 3
             if (nside2npix(nside_in) == nsize) ftype_in = 2
             goto 1000
          endif
       endif
1000   continue

       ! find out if polarisation data is present
       if (present(polarisation)) then
          if (tfields < 3) then
             polarisation = 0 ! no polarisation
             goto 2000
          endif
          call ftgkyl(unit,'POLAR',polar_in,comment,status)
          if (status == 0) then ! polar found
             polarisation = 0
             if (polar_in) polarisation = 1
             goto 2000
          else ! polar not found
             status = 0
             polarisation = -1
             if (hdutype <= 0) goto 2000
             if (hdutype == 1) then ! ascii table -> power spectra
                ndp = 4
                defpol(1:ndp,1) = (/ "GRA","E-M","POW","EE " /)
                defpol(1:ndp,2) = (/ "CUR","B-M","POW","BB " /)    
             endif
             if (hdutype == 2) then ! binary table -> maps or power spectra
                ndp = 10
                defpol(1:ndp,1) = (/ "Q-P","Q_P","Q P", "Q-S","Q_S", &
                     &               "Q  ",      "GRA","E-M","POW","EE " /)
                defpol(1:ndp,2) = (/ "U-P","U_P","U P", "U-S","U_S", &
                     &               "U  ",      "CUR","B-M","POW","BB " /)
             endif
             pf(:) = .false.
             do i = 2, tfields ! do not consider first field (generally temperature)
                ttype_val = adjustl(ttype(i))
                call ftupch(ttype_val) ! upper case
                do k=1,2
                   do j=1, ndp
                      if (index( ttype_val, trim(defpol(j,k)) ) == 1) then
                         pf(k) = .true.
                         goto 1500 ! move to next field
                      endif
                   enddo
                enddo ! loop on k
1500            continue
             enddo ! loop on i
             if (pf(1) .and. pf(2)) polarisation = 1
             if (.not. (pf(1) .or. pf(2))) polarisation = 0
          endif ! polar not found
       endif ! present(polarisation)
2000   continue

       call ftgkys(unit,'POLCCONV',polcconv_val,comment,status)
       if (status == 0) then
          polcconv_in = 3 ! neither COSMO nor IAU
          if (trim(polcconv_val) == 'COSMO') polcconv_in = 1
          if (trim(polcconv_val) == 'IAU')   polcconv_in = 2
       else
          polcconv_in = 0 ! not found
          status = 0
       endif

    else ! no image no extension, you are dead, man
       ftype_in = -1
       call fatal_error(' No image, no extension')
    endif

    !     close the file
    call ftclos(unit, status)

    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)

    getsize_fits=nsize

    call ftupch(order_val) ! convert order_val to upper case
    if (order_val(1:4) == 'RING') ordering_in = 1
    if (order_val(1:4) == 'NEST') ordering_in = 2

    if (present(nmaps))       nmaps       = nmaps_in
    if (present(mlpol))       mlpol       = mlpol_in
    if (present(obs_npix))    obs_npix    = obs_npix_in
    if (present(ordering))    ordering    = ordering_in
    if (present(nside))       nside       = nside_in
    if (present(type))        type        = ftype_in
    if (present(fwhm_arcmin)) fwhm_arcmin = fwhm_arcmin_in
    if (present(beam_leg))    beam_leg    = adjustl(beam_leg_in)
    if (present(coordsys))    coordsys    = adjustl(coordsys_in)
    if (present(polcconv))    polcconv    = polcconv_in

    return
  end function getsize_fits
  !=======================================================================
  function number_of_alms(filename, extnum)
    !=======================================================================
    !    Read the number of alms from a FITS-file containing alms          
    !    for constrained realisations in synfast.
    !
    !         Frode K. Hansen, April 1999
    !    EH. Jan 2004: use repeat information
    !=======================================================================
    character(len=*),           intent(in)    :: filename
    integer(I4B),     optional, intent(inout) :: extnum
    integer(I4B)                              :: number_of_alms
    !integer(I8B)                              :: number_of_alms

    integer(I4B), dimension(2)           :: naxes
    integer(I4B) :: status, unit, readwrite, blocksize, naxis, nfound
    integer(I4B) :: nmove, hdutype, hdunum
    integer(I4B) :: datacode, repeat, width
    character(len=80) :: comment
    character(len=20) :: tform
    logical(LGT)      ::  extend

    !-----------------------------------------------------------------------
    status=0
    unit = 118
    comment=''

    readwrite=0
    call ftopen(unit, filename, readwrite, blocksize, status)
    if (status > 0) call printerror(status)
    !     -----------------------------------------
    call ftgkyj(unit,'NAXIS', naxis, comment, status)

    !     determines the presence of an extension
    call ftgkyl(unit,'EXTEND', extend, comment, status)
    call assert (status<=0, 'No Extension in FITS file!')

    nmove = +1
    call ftmrhd(unit, nmove, hdutype, status)

    call assert (hdutype==2, 'This is not a FITS binary-table')

    call ftgknj(unit,'NAXIS', 1_i4b, 2_i4b, naxes, nfound, status)
    call assert (nfound>=2, 'NAXIS2-keyword not found!')

    call ftgkys(unit,'TFORM1', tform, comment, status)
    call ftbnfm(tform, datacode, repeat, width, status)

    number_of_alms = int(naxes(2),kind=i8b) * repeat

    if (present(extnum)) then
       call ftthdu(unit,hdunum,status)
       extnum = hdunum - 1 ! first HDU : primary array
    endif

    call ftclos(unit, status)

  end function number_of_alms

  !======================================================================
  subroutine putrec(unit, card, status)
    !======================================================================
    ! append, delete or update a card from the current header
    ! (deletion of keyword KWD is described by : - KWD)
    ! EH, version 1.0, Dec 2001
    !======================================================================
    integer(kind=I4B), intent(IN)  :: unit
    character(len=*),  intent(IN)  :: card
    integer(kind=I4B), intent(OUT) :: status

    character(len=80) :: cardfits,record
    character(len=8)  :: kwd
    integer(kind=I4B) :: hdtype
    character(len=80) :: nullarr(0)
    !======================================================================

    status = 0
    cardfits=''
    record=''
    call ftgthd(card, cardfits, hdtype, status)
    kwd = cardfits(1:8)
    status = 0

    select case (hdtype)
    case (-1) ! delete keyword (starting from the top)
       call ftgrec(unit,0_i4b,record,status)
       !             call ftdkey(unit, kwd, status)
       ! patch around cfitsio bug 
       ! (ftdkey does not recognize wild cards while ftgnxk does)
       do
          call ftgnxk(unit, kwd, 1_i4b, nullarr, 0_i4b, record, status)
          if (status /= 0) exit
          call ftdkey(unit, record(1:8), status)
       enddo
    case (0) ! append or update
       if (kwd == 'CONTINUE') then
          call ftprec(unit, trim(card), status)
       else
          ! if long string, put LongStringWarning
          if (index(card, "&'")>0) call ftplsw(unit, status)
          ! delete keyword in its current location (if any)
          call ftdkey(unit, kwd, status)
          status = 0
          ! append
          call ftprec(unit, trim(cardfits), status)
       endif
    case (1) ! append (for HISTORY and COMMENT)
       call ftprec(unit, trim(cardfits), status)
    case default
       write(unit=*,fmt=*)" Unexpected card format in fits header :"
       write(unit=*,fmt="(a80)") card
       write(unit=*,fmt=*)" card not written."
    end select
    status = 0

    return
  end subroutine putrec

  !====================================================================
  subroutine get_clean_header(unit, header, filename, error, xalso, xonly)
    !====================================================================
    ! get_clean_header(unit, header, error [, xalso, xonly])
    ! gets the FITS header from unit, after excluding some default keywords 
    !  defined in def_excl
    ! if header in non void on input, its content will be concatenated with that
    !  of the FITS header
    ! if xalso is defined as a list of keywords, they are also excluded from the header
    ! if xonly is defined as a list of keywords, only those keywords are excluded from
    ! the header.
    ! xonly and xalso are exclusive
    !====================================================================
    INTEGER(I4B),                    intent(IN)           :: unit
    CHARACTER(LEN=*), DIMENSION(1:), INTENT(IN OUT)       :: header
    CHARACTER(LEN=*),                INTENT(IN)           :: filename
    INTEGER(I4B),                    intent(OUT)          :: error
    character(len=8), dimension(1:), intent(IN), optional :: xalso
    character(len=8), dimension(1:), intent(IN), optional :: xonly

    INTEGER(I4B) :: nlheader, status, i, n_excl
    CHARACTER(LEN=80) :: record
    CHARACTER(len=8), dimension(:), allocatable :: to_excl

    CHARACTER(len=8), dimension(1:21) :: def_excl
    !====================================================================

    ! keywords to be excluded by default from output header
    ! Note that TTYPE# and TUNIT# keywords are not excluded
    ! from the output header as they might be useful to the 
    ! calling routines
    def_excl=(/&
         & "SIMPLE  ","BITPIX  ","NAXIS   ",&
         & "NAXIS#  ","PCOUNT  ","GCOUNT  ",&
         & "EXTEND  ","ORIGIN  ","DATE*   ",&
         & "TFIELDS ","TFORM#  ",           & 
         & "TBCOL#  ","EXTNAME ","CTYPE#  ",&
         & "CRVAL#  ","CRPIX#  ","CDELT#  ",&
         & "XTENSION","INSTRUME","TELESCOP",&
         & "PDMTYPE "/)

    error = 0
    record=''

    if (present(xonly)) then 
       n_excl = size(xonly)
       allocate(to_excl(1:n_excl))
       to_excl = xonly

    else if (present(xalso)) then
       n_excl = size(xalso) + size(def_excl)
       allocate(to_excl(1:n_excl))
       to_excl(1:size(def_excl)) = def_excl
       to_excl(size(def_excl)+1:n_excl) = xalso

    else
       n_excl = size(def_excl)
       allocate(to_excl(1:n_excl))
       to_excl = def_excl
    endif

    nlheader=size(header)
    ! go to end of fortran header
    do i = 1, nlheader
       if (trim(header(i)) == "") exit
    enddo
    ! go to top of fits file header
    status=0
    call ftgrec(unit,0_i4b,record,status)
    ! read in all header lines except those excluded
    do
       call ftgnxk(unit,'*',1_i4b,to_excl,n_excl,record,status)
       if (status > 0) exit ! end of header
       if (i > nlheader) then
          write(unit=*,fmt="(a,i5,a)") &
               & " WARNING : The header in "//  &
               &    trim(filename)//" has more than ", &
               &  nlheader," lines."
          print*," It will be truncated."
          error = 1
          exit
       endif
       header(i)=record
       i=i+1
    enddo
    status=0

    return
  end subroutine get_clean_header

  
  !====================================================================
  subroutine check_input_map(code, mapfile, polarisation)
  !====================================================================
    ! check out that file contains a valid HEALPIX map
    ! on input polarisation is set to T if one wants to get polarisation data from map,
    ! and on output it will be set to F if that is not possible
    !====================================================================
    use pix_tools, only: nside2npix, npix2nside
    character(len=*)                      :: code
    character(len=*)                      :: mapfile
    logical(LGT), intent(inout)           :: polarisation
    
    !
    integer(I4B) :: nmaps, order_map, nsmax, mlpol, type, polar_fits, polcconv
    integer(I4B) :: npixtot
    character(len=1) :: coordsys
    character(len=*), parameter :: primer_url = 'http://healpix.sf.net/pdf/intro.pdf'
    !====================================================================
    
    !n_ext = getnumext_fits(mapfile)
    npixtot = getsize_fits(mapfile, nmaps = nmaps, ordering=order_map, nside=nsmax,&    
         &              mlpol=mlpol, type = type, polarisation = polar_fits, &
         &             coordsys=coordsys, polcconv=polcconv)
    
    if (nsmax<=0) then
       print*,"Keyword NSIDE not found in FITS header of "//trim(mapfile)//" !"
       call fatal_error(code)
    endif
    if (type == 3) npixtot = nside2npix(nsmax) ! cut sky input data set
    if (nsmax/=npix2nside(npixtot)) then
       print 9000,"FITS header keyword NSIDE does not correspond"
       print 9000,"to the size of the map in "//trim(mapfile)//" !"
       call fatal_error(code)
    endif

    if (polarisation .and. (nmaps >=3) .and. polar_fits == -1) then
       print 9000,"The input FITS file "//trim(mapfile)//" MAY NOT contain polarisation data."
       print 9000,"Proceed at your own risk"
    endif
    
    if (polarisation .and. (nmaps<3 .or. polar_fits ==0)) then
       print 9000,"The FITS file "//trim(mapfile)//" does NOT contain polarisation maps"
       print 9000,"only the temperature field will be analyzed"
       polarisation = .false.
    endif
    
    ! POLCCONV now dealt with in input_map
!     if (polarisation .and. (polcconv == 3)) then
!        print 9000,"The polarisation coordinate convention (POLCCONV) is neither COSMO nor IAU."
!        print 9000,code//" can not proceed with these data"
!        print 9000,"See Healpix primer ("//primer_url//") for details."
!        call fatal_error(code)
!     endif
    
!     if (polarisation .and. (polcconv == 2)) then
!        print 9000,"The input map contains polarized data in the IAU coordinate convention (POLCCONV)"
!        print 9000,code//" can not proceed with these data"
!        print 9000,"See Healpix primer ("//primer_url//") for details."
!        call fatal_error(code)
!     endif
    
!     if (polarisation .and. (polcconv == 0)) then
!        print 9000,"WARNING: the polarisation coordinate convention (POLCCONV) can not be determined"
!        print 9000,"         COSMO will be assumed."
!        print 9000,"See Healpix primer ("//primer_url//") for details."
!     endif
    
    !     --- check ordering scheme ---
    if ((order_map/=1).and.(order_map/=2)) then
       print 9000,"The ordering scheme of the map must be RING or NESTED."
       print 9000,"No ordering specification is given in the FITS-header of "//trim(mapfile)//" !"
       call fatal_error(code)
    endif
    
9000 format(a)
    
    return
  end subroutine check_input_map
      

end module fitstools
