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

module mask_tools
  ! module for mask processing
  !
  ! All input mask are assumed to be in NESTED ordering
  !
  ! HISTORY:
  !  E. Hivon 2006-2010 firt developments (mask_processing.f90)
  !  2011-11-22: inclusion in Healpix 

  use healpix_types
  use long_intrinsic, only: long_size
  implicit none

  private

  public :: fill_holes_nest
  public :: maskborder_nest
  public :: dist2holes_nest
  public :: size_holes_nest

  integer(i4b), parameter, private :: PIXBAD = 0  ! invalid pixels
  integer(i4b), parameter, private :: PIXGOOD = 1 ! valid pixels
  integer(i4b), parameter, private :: DEFBTAG = 2 ! default border value, must not be 0 nor 1


contains
  !======================================================================
  subroutine fill_holes_nest(nside, newminsize, mask_in, mask_out)
    !======================================================================
    ! fill_holes_nest(nside, newminsize, mask_in, mask_out)
    !
    ! nside      IN, Healpix resolution parameter
    ! newminsize IN, minimal size of hole (in pixels) on output
    ! mask_in    IN, Integer*4, full sky map, in NESTED scheme
    !       mask of good pixels (mask=1) and bad pixels (mask=0)
    ! mask_out   OUT, Integer*4, full sky map, in NESTED scheme
    !       bad pixels located in a 'hole' with fewer than newminsize pixels
    !      take value 1 (ie become valid)
    !
    ! A hole is the set of all adjacent pixels initially set to 0.
    ! 2 pixels are adjacent if they have at least one point in common.
    !
    !======================================================================
    use pix_tools,  only: nside2npix
    use misc_utils, only: fatal_error
    ! dummy variables
    integer(I4B),                intent(IN)  :: nside
    integer(I4B),                intent(IN)  :: newminsize
    integer(I4B), dimension(0:), intent(IN)  :: mask_in
    integer(I4B), dimension(0:), intent(OUT) :: mask_out
    ! local variables
    integer(I4B) :: out
    integer(I4B), parameter :: MKD = I4B
    integer(MKD) :: p, nholes, npholes
    integer(MKD), dimension(:), pointer :: sizeholes, listpix
    !-------------------------------------------------------------------

    mask_out = mask_in
    if (newminsize > 0) then
       call size_holes_nest(nside, mask_out, nholes, npholes, sizeholes=sizeholes, listpix=listpix)
       do p=0,nholes-1
          out = PIXBAD ! leave that hole to hole value
          if (sizeholes(p) < newminsize) out = PIXGOOD ! fill that hole
          mask_out(listpix(listpix(p):listpix(p+1)-1)) = out
       enddo
    endif

    return
  end subroutine fill_holes_nest

  !======================================================================
  subroutine maskborder_nest(nside, mask_in, mask_out, nbordpix, border_value)
  !======================================================================
    ! maskborder_nest(nside, mask_in, mask_out, nbordpix [, border_value])
    !
    ! mask_in IN, Integer*4, full sky map, in NESTED scheme
    !    on input: mask of good pixels (mask=1) and bad pixels (mask=0)

    ! mask_out OUT, Integer*4, full sky map, in NESTED scheme
    !    on output: inner border of bad regions takesmask=border_value
    !
    ! nbordpix OUT, Integer*4
    !    number of pixels in border
    !
    ! border_value IN, OPTIONAL, Integer*4
    !    value given to border pixels = 2 by default
    !======================================================================
    use pix_tools,  only: nside2npix, neighbours_nest
    use misc_utils, only: fatal_error

    ! dummy variables
    integer(I4B),                intent(IN) :: nside
    integer(I4B), dimension(0:), intent(in) :: mask_in
    integer(I4B), dimension(0:), intent(out) :: mask_out
    integer(I4B),                intent(out) :: nbordpix
    integer(I4B), optional,      intent(in) :: border_value
    
    ! local variables
    integer(I8B) :: npix, sz1, sz2, p
    integer(I8B), dimension(0:7) :: list
    integer(I4B) :: minmask, maxmask, i, nn
    integer(I4B), dimension(0:7) :: values
    character(len=*), parameter  :: code = 'maskborder_nest'
    integer(i4b) :: mytag ! user-defined border value
    integer(i4b) :: minallowed, maxallowed
    !----------------------------------------------------------------------

    sz1  = long_size(mask_in)
    sz2  = long_size(mask_out)
    npix = nside2npix(nside)
    if (sz1 /= npix .or. sz2 /= npix) then
       write(*,'(a,a,i11,a)')    code,' found ',   sz1,  ' pixels in mask_in'
       write(*,'(a,a,i11,a)')    code,' found ',   sz2,  ' pixels in mask_out'
       write(*,'(a,a,i11,a,i5)') code,' expected ',npix,' pixels for Nside = ',nside
       call fatal_error('Aborting '//code)
    endif
    minmask = minval(mask_in)
    maxmask = maxval(mask_in)
    minallowed = min(PIXGOOD, PIXBAD)
    maxallowed = max(PIXGOOD, PIXBAD)
    if (maxmask > maxallowed .or. minmask < minallowed) then
       print*,'>>>>> Range of input mask: ',minmask,maxmask
       print*,'>>>>> Allowed range: ',minallowed, maxallowed
       call fatal_error('Aborting '//code)
    endif

    ! test particular cases with no border
    nbordpix = 0
    if (maxmask == minallowed) return !all pixels are bad
    if (minmask == maxallowed) return !all pixels are good

    mask_out = mask_in
    ! start actual calculation
    do p = 0, npix - 1
       if (mask_out(p) == PIXBAD) then
          ! test neighbours of bad pixels
          call neighbours_nest(nside, p, list, nn)
          values(0:nn-1) = mask_out(list(0:nn-1))
          do i=0,nn-1
             ! if any neighbour is valid, then bad pixel belongs to border
             if (values(i) == PIXGOOD) then
                nbordpix = nbordpix + 1
                mask_out(p) = DEFBTAG
                goto 1234
             endif
          enddo
1234      continue

       endif
    enddo

    ! change border value from default to user-defined (if applicable)
    mytag = DEFBTAG
    if (present(border_value)) then
       mytag = border_value
    endif
    if (mytag /= DEFBTAG) then
       do p = 0, npix - 1
          if (mask_out(p) == DEFBTAG) mask_out(p) = mytag
       enddo
    endif

    return
  end subroutine maskborder_nest

  ! ==========================================================================
  subroutine dist2holes_nest(nside, mask, distance)
    ! ==========================================================================
    ! for each valid pixel, returns distance (in radians) from pixel center to 
    ! center of closest invalid pixel.
    ! The input mask must be NESTED ordered
    !
    ! Algorithm:
    !  To begin with, the inner border (=invalid pixels in touch with valid ones) 
    !  of each invalid region is identified
    ! Then a small Nside is used to select the closest border pixels, and finally
    ! the distance of each valid pixel to the selected border pixels is computed
    ! to find its minimum.
    !
    ! Note: final distances are computed with  2 ASIN (|v1-v2|/2) which is more
    ! accurate than the usual ACOS(v1.v2) at small separations.
    ! ==========================================================================

    use pix_tools, only: nside2npix, npix2nside, pix2vec_nest
    use misc_utils, only: fatal_error, wall_clock_time
    use num_rec, only: indexx

    integer(I4B),                intent(IN)  :: nside
    integer(I4B), dimension(0:), intent(IN)  :: mask
    real(DP),     dimension(0:), intent(OUT) :: distance
!    integer(i4b), parameter :: algo = 1
    integer(i4b), parameter :: nslow = 64 ! 64 or 128
!    integer(i4b), parameter :: dimtest = 1
    logical(LGT), parameter :: brute_force = .false.
    logical(LGT), parameter :: do_clock = .false. !.true.
    

    integer(I8B) :: npix1, npix2, npix, p, iratio
    integer(I4B) :: nside1, nside2
    integer(I4B) :: pb, pv, q, qq
    integer(I4B) :: cache, nstride, nl, il, i1, i2
    integer(I4B) :: n_invalid
    integer(I4B), dimension(:), allocatable   :: mask1
    integer(I4B)                              :: nbordpix
    real(DP),     dimension(:,:), allocatable :: bordpos
    real(DP),     dimension(:,:), allocatable :: bv
    real(DP),     dimension(:  ), allocatable :: dd
    integer(I4B), dimension(:  ), allocatable :: bordpix
    real(DP),     dimension(1:3)              :: vv, v1
    real(DP)                                  :: ddmin
    integer(I4B) :: nslow_in, nplow, nploweff, n2test, n_in_list
    integer(I4B) :: nlactive
    real(DP),    dimension(:,:,:), allocatable :: vertlow
    integer(I4B), dimension(:), allocatable :: nb_in_lp, p2test, start_list
    integer(I4B), dimension(:,:), allocatable :: non_empty_low
    real(DP) :: ddmax, dtotop
    real(DP), dimension(1:3,1:4) :: tmpvertex

    real(DP), dimension(:), allocatable :: drange
    real(DP)                              :: distmin, distmax

    real(DP),    dimension(:,:), allocatable :: centlow, centpix

    character(len=*), parameter :: code = 'dist2holes_nest'
    integer(i4b) ::  algo_in, q1, pp, pstart, np, nb
    real(DP),     dimension(:,:), allocatable :: bordpospack
    real(sp) :: t0, t1, t2, t3, ts1, ts2, ts3, tstart, tend, ta1, tsa1, tinit
    real(DP) :: fudge, cfudge, sfudge, threshold, mean_nbord
    !----------------------------------------------------------------------

    ! test inputs
    npix1 = long_size(mask)
    nside1 = npix2nside(npix1)
    if (nside1 < 0) call fatal_error('Invalid size of input mask in '//code)

    npix2 = long_size(distance)
    nside2 = npix2nside(npix2)
    if (nside2 < 0) call fatal_error('Invalid size of distance in '//code)

    npix = nside2npix(nside)
    if (npix1 /= npix .or. npix2 /= npix) then
       call fatal_error('Arrays too small in '//code)
    endif

    ! find pixels lying on holes border
    allocate(mask1(0:npix-1))
    call maskborder_nest(nside, mask, mask1, nbordpix)

    ! treat particular case of no border
    if (nbordpix == 0) then
       distance = 0.0_DP
       return
    endif
    n_invalid = 0

    !algo_in = algo

    if (do_clock) call wall_clock_time(tstart)
    pb = 0
    distance = 1.e6 ! set all pixels to large values
    ! find 3D location of those border pixels
    allocate(bordpos(0:nbordpix-1,1:3))
    allocate(bordpix(0:nbordpix-1))
    do p=0, npix-1
       if (mask1(p) == DEFBTAG) then 
          call pix2vec_nest(nside, p, vv)
          bordpos(pb,1:3) = vv
          bordpix(pb) = p
          pb = pb + 1
          distance(p) = 0.0_dp
       endif
       if (mask1(p) == PIXBAD) distance(p) = 0.0_dp
    enddo

    !if (algo_in == 1) then
       ! define low resolution grid
       nslow_in = min(nslow, nside)
       nplow = nside2npix(nslow_in)
       iratio = nside / nslow_in
       iratio = iratio * iratio ! ratio of number of pixels
       fudge = 4.276_dp / real(nslow_in, kind=DP) ! 2.0 x diagonal of longest pixel, at large Nside
       fudge = min(fudge, PI)
       sfudge = sin(fudge)
       cfudge = cos(fudge) ! cosine of fudge

       ! find low resolution pixels containing border pixels (=non_empty)
       allocate(nb_in_lp(0:nplow-1))
       nb_in_lp = 0
       do p=0, nbordpix-1
          q = bordpix(p) / iratio
          nb_in_lp(q) = nb_in_lp(q) + 1 ! number of border pixels in each low-res pixel
       enddo
       if (sum(nb_in_lp) /= nbordpix) then
          print*,sum(nb_in_lp),nbordpix
          call fatal_error('Error 1 with low res pixels')
       endif
       if (maxval(nb_in_lp) > iratio) then
          print*,maxval(nb_in_lp),iratio
          call fatal_error('Error 1a with low res pixels')
       endif
       nploweff = count(nb_in_lp > 0)
       allocate(non_empty_low(0:nploweff-1,1:2))
       qq = 0
       do q=0, nplow-1
          if (nb_in_lp(q) > 0) then
             non_empty_low(qq,1) = q ! list of non-empty low-res pixels
             non_empty_low(qq,2) = nb_in_lp(q) ! number of border pixels
             qq = qq + 1
          endif
       enddo
       allocate(start_list(0:nploweff))
       start_list(0) = 0
       do qq=1, nploweff
          start_list(qq) = start_list(qq-1) + non_empty_low(qq-1,2)
       enddo
       if (start_list(nploweff) /= nbordpix) then
          print*,start_list(nploweff),nbordpix
          call fatal_error('Error 2 with low res pixels')
       endif
       deallocate(nb_in_lp)

       ! build location of center for non-empty low resolution "border" pixels
       allocate(centlow(1:3,    0:nploweff-1))
       do qq=0, nploweff-1
          q = non_empty_low(qq,1)
          call pix2vec_nest(nslow_in, q, vv)
          centlow(1:3, qq)      = vv(1:3)
       enddo
       allocate(drange(0:nploweff-1))
       allocate(p2test(    0:nploweff-1))
       allocate(centpix(1:3, 0:iratio-1))
       allocate(dd(0:iratio-1))
       allocate(bordpospack(0:nbordpix-1, 1:3 ))

       if (do_clock) then
          call wall_clock_time(t0)
          tinit = t0-tstart
       endif

       ts1 = 0. ; tsa1 = 0. ; ts2 = 0. ; ts3 = 0.
       nlactive = 0
       ! loop on low-resolution pixels
       do q1=0, nplow-1
          ! keep only those containing valid small pixel
          if (any(mask1(q1*iratio:(q1+1)*iratio-1) == PIXGOOD) ) then
             nlactive = nlactive + 1
             if (do_clock) call wall_clock_time(t0)
             ! center location of current low-res pixel
             call pix2vec_nest(nslow_in, q1, vv)

             ! distance of low-res pixel center with each center of low-res "border" pixels
! !$OMP PARALLEL DEFAULT(NONE) &
! !$OMP SHARED(nploweff, drange, centlow, vv) &
! !$OMP PRIVATE(qq)
! !$OMP DO schedule(static, 4)
             do qq=0, nploweff - 1
                drange(qq) = centlow(1, qq) * vv(1) &
                     &     + centlow(2, qq) * vv(2) &
                     &     + centlow(3, qq) * vv(3) 
             enddo
! !$OMP END DO
! !$OMP END PARALLEL
             ! min distance (=max scalar product)
             distmin = maxval(drange)
             distmin = min(distmin, 1.0_dp)
             
             ! find out list of low-res "border" pixels close enough to be checked-out at higher-res
             ! radius = smallest distance found above + fudge
             if ((distmin + cfudge) > 0.0_dp) then ! new radius is less than Pi
                threshold = distmin*cfudge - sfudge*sqrt((1.d0-distmin)*(1.d0+distmin))
             else! new radius can not exceed Pi
                threshold = -1.0_dp
             endif
             n2test = 0
             do qq=0, nploweff - 1
                if (drange(qq) >= threshold) then ! for scalar product
                   p2test(n2test) = qq
                   n2test = n2test + 1
                endif
             enddo
                
             ! store location of border pixels contained in relevant low-res pixels
             pstart = 0
             do q=0, n2test-1
                qq = p2test(q)
                i1 = start_list(qq)
                i2 = start_list(qq+1)-1
                n_in_list = i2 - i1 + 1
                bordpospack(pstart:pstart+n_in_list-1,1:3) = bordpos(i1:i2,1:3)
                pstart = pstart + n_in_list
             enddo
             nb = pstart ! number of border pixels close enough to be considered
             mean_nbord = mean_nbord + nb

             if (do_clock) then
                call wall_clock_time(t1)
                ts1 = ts1 + (t1 - t0)
             endif
             ! center location for small pixels
             pp = 0
             do p=q1*iratio, (q1+1)*iratio - 1
                if (mask1(p) == PIXGOOD) then
                   call pix2vec_nest(nside, p, v1)
                   centpix(1:3, pp) = v1
                   pp = pp + 1
                endif
             enddo
             np = pp ! number of valid small pixels in current low-res pixels
             
             if (do_clock) then
                call wall_clock_time(ta1)
                tsa1 = tsa1 + (ta1 - t1)
             endif
             ! scalar product: small pixels * border pixels

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(np, nb, dd, bordpospack, centpix) &
!$OMP PRIVATE(pp)
!$OMP DO schedule(guided, 4)
             do pp = 0, np-1
                ! compute (v1-v2)^2 in [0,4]
                dd(pp) = minval(  (bordpospack(0:nb-1,1)-centpix(1,pp))**2 &
                     &          + (bordpospack(0:nb-1,2)-centpix(2,pp))**2 &
                     &          + (bordpospack(0:nb-1,3)-centpix(3,pp))**2)
             enddo
!$OMP END DO
!$OMP END PARALLEL

             if (do_clock) then
                call wall_clock_time(t2)
                ts2 = ts2 + (t2 - ta1)
             endif

             ! norm2 -> angular distance (radians)
             pp = 0
             do p=q1*iratio, (q1+1)*iratio - 1
                if (mask1(p) == PIXGOOD) then
                   if (dd(pp) < 0.0_dp .or. dd(pp) > 4.0_dp) then
                      n_invalid = n_invalid + 1
                   endif
                   distance(p) = 2.d0 * asin( sqrt(dd(pp)) * 0.5d0 )
                   pp = pp + 1
                endif
             enddo
             if (do_clock) then
                call wall_clock_time(t3)
                ts3 = ts3 + (t3 - t2)
             endif
          endif
       enddo ! loop on low res pixels
       deallocate(drange, p2test, dd)
       deallocate(bordpix)
       deallocate(non_empty_low)
       if (n_invalid > 0) then
          write(*,'(a,i7,a)') 'WARNING: Found ',n_invalid,' invalid distances: '
       endif
       if (do_clock) then
          call wall_clock_time(tend)
          print*,nside, nslow, nbordpix, mean_nbord/nlactive
          write(*,'(7(a12))')     'init', 'prep', 'pix2vec', 'product','acos', 'total',          'actual'
          write(*,'(4x,7(g12.4))') tinit, ts1,     tsa1,      ts2,    ts3,    ts1+tsa1+ts2+ts3, tend-tstart
       endif

    !endif ! algo_in
       
    deallocate(mask1)

  end subroutine dist2holes_nest

  !=========================================================================
   subroutine size_holes_nest(nside, mask, nholes, nph, &
        &                     tags, sizeholes, listpix)
     !=========================================================================
     ! SIZE_HOLES_NEST: subroutine to determine size (=number of pixels) of holes
     ! A hole is the set of all adjacent pixels initially set to 0.
     ! 2 pixels are adjacent if they have at least one point in common.
     !
     ! Nside: integer, IN, resolution parameter
     ! Mask(0:): integer 1D array, IN, each pixel must be either 1 (=valid) or 0 (=invalid)
     ! Nholes: integer, OUT, number of holes
     ! Nph:    integer, OUT, total number of pixels in holes
     ! Tags(0:npix-1): integer 1D array, OPTIONAL, OUT:
     !                   invalid pixels belonging to largest hole have value -1,
     !                   invalid pixels belonging to second largest hole: -2,
     !                   and so on, while valid pixels have value 1
     ! Sizeholes(0:Nholes-1): integer pointer, OPTIONAL, OUT: respective size of each hole
     ! Listpix(0:Nph+Nholes): integer pointer, OPTIONAL, OUT: list of pixels in each hole
     !
     !========================================================================
    use healpix_types
    use pix_tools,  only: nside2npix, neighbours_nest
    use misc_utils, only: fatal_error, assert_alloc
    use num_rec,    only: iindexx

    integer(I4B),                intent(IN) :: nside
    integer(I4B), dimension(0:), intent(IN) :: mask
    integer(I4B),                intent(OUT)   :: nholes, nph
    integer(I4B), dimension(:),  optional, pointer :: tags
    integer(I4B), dimension(:),  optional, pointer :: sizeholes
    integer(I4B), dimension(:),  optional, pointer :: listpix

    integer(I4B) :: i, p, q, nnb, nvalid, nbad, status
    integer(I4B) :: minmask, maxmask, nlist
    integer(I4B) :: h, k, k2, np, shift
    integer(I4B), dimension(0:7) :: lnb
    integer(I4B) :: ninstack, stackmax, hole_id
    integer(I4B), dimension(:), allocatable :: stack, list, mytags
    integer(I4B), dimension(:), allocatable :: tally
    integer(I4B), dimension(:), allocatable :: indx, cumulu, cumuls, revert
    integer(I4B), parameter :: unprocessed = 0
    integer(I8B) :: npix, sz
    character(len=*), parameter  :: code = 'size_holes_nest'
    logical(LGT), parameter :: sorted = .true.
    logical(LGT) :: do_size, do_list, do_tags
    !----------------------------------------------------------------------

    sz = size(mask,1)
    npix = nside2npix(nside)
    if (sz /= npix) then
       write(*,'(a,a,i11,a)')    code,' found ',   sz,  ' pixels in mask'
       write(*,'(a,a,i11,a,i5)') code,' expected ',npix,' pixels for Nside = ',nside
       call fatal_error('Aborting '//code)
    endif
    minmask = minval(mask)
    maxmask = maxval(mask)
    if (maxmask > 1 .or. minmask < 0) then
       print*,'>>>>> Range: ',minmask,maxmask
       call fatal_error('input mask should be in [0,1] in '//code)
    endif

    do_size = (present(sizeholes))
    do_list = (present(listpix))
    do_tags = (present(tags))

    nvalid = sum(mask)
    nbad = npix - nvalid
    stackmax = max(8*nbad,100)

    allocate(stack(0:stackmax), stat=status)
    call assert_alloc(status, code, 'stack')
    allocate(tally(0:npix/2), stat=status)
    call assert_alloc(status, code, 'tally')
    allocate(list(0:nbad+1), stat=status)
    call assert_alloc(status, code, 'list')
    allocate(mytags(0:npix-1), stat=status)
    call assert_alloc(status, code, 'mytags')
    stack(:) = -1
    tally(:) = 0
    list(:)  = -1
    k        = 0
    hole_id  = 0
    mytags   = mask

    !print*,npix, nvalid, nbad


    ! loop over all pixels
    do p=0, npix-1

       if (mytags(p) == unprocessed) then
          ! start a new hole and a new stack
          hole_id = hole_id + 1
          ninstack = 1
          q = p

          do

             ! mark current pixel as processed
             ! add it to current hole tally and hole list
             if (mytags(q) == unprocessed) then
                mytags(q) = -hole_id
                tally(hole_id-1) = tally(hole_id-1) + 1
                list(k) = q
                k = k + 1
             endif

             ! remove it from stack
             ninstack = ninstack - 1

             ! put unprocessed neighbouring pixels 
             !    belonging to hole into stack
             call neighbours_nest(nside, q, lnb, nnb)
             do i=0,nnb-1
                if (mytags(lnb(i)) == unprocessed) then
                   stack(ninstack) = lnb(i)
                   ninstack = ninstack + 1
                endif
             enddo

             ! move on to pixel on top of stack, and repeat
             !    unless stack is empty
             if (ninstack == 0) exit
             q = stack(ninstack-1)
          enddo

       endif

    enddo
    deallocate(stack)

    ! create OPTIONAL Tags output
    if (do_tags) then 
       allocate(tags(0:npix-1),stat=status)
       call assert_alloc(status, code, 'tags')
       tags(0:npix-1) = mytags(0:npix-1)
    endif
    deallocate(mytags)


!     print*,'number of holes:',hole_id
!     if (hole_id > 0) then
!        print*,'Hole,    size'
!        do p=0, hole_id-1
!           print*,p,tally(p)
!        enddo
!     endif
!     print*,'total',nbad,sum(tally)

    nholes = hole_id
    nph = sum(tally)
    nlist = nholes + nph + 1


    if (do_size .or. do_list .or. (do_tags.and.sorted)) then
       if (.not. sorted) then
          ! ----- unsorted lists --------
          if (do_size) then
             allocate(sizeholes(0:nholes-1), stat=status)
             call assert_alloc(status, code, 'sizeholes')
             sizeholes = tally
          endif
          
          if (do_list) then
             allocate(listpix(0:nlist-1), stat=status)
             call assert_alloc(status, code, 'listpix')

             allocate(cumulu(0:nholes), stat=status)
             call assert_alloc(status, code, 'cumulu')
             cumulu(0) = 0
             do h=1,nholes
                cumulu(h) = cumulu(h-1) + tally(h-1)
             enddo

             listpix(0:nholes) = cumulu+nholes+1
             listpix(nholes+1:nlist-1) = list
             deallocate(cumulu)
          endif

       else
          ! ----- sorted lists -------

          allocate(indx(0:nholes-1), stat=status)
          call assert_alloc(status, code, 'indx')
          ! find indexing for size sorting (decreasing order)
          tally = -tally
          call iindexx(nholes, tally, indx)
          tally = -tally
          shift = minval(indx) - lbound(tally,1)
          indx = indx - shift ! make sure INDX is consistent with TALLY indexing

          if (do_size) then
             allocate(sizeholes(0:nholes-1), stat=status)
             call assert_alloc(status, code, 'sizeholes')
             sizeholes = tally(indx)
          endif

          if (do_tags) then
             allocate(revert(0:nholes-1), stat=status)
             call assert_alloc(status, code, 'revert')
             do h=0, nholes-1
                revert(indx(h)) = h ! in 0,1,2, ..
             enddo
             do p=0,npix-1
                q = tags(p) ! in -1, -2, -3, ..
                if (q < 0) then
                   tags(p) = -1 - revert(-1-q)
                endif
             enddo
             deallocate(revert)
          endif

          if (do_list) then
             allocate(listpix(0:nlist-1), stat=status)
             call assert_alloc(status, code, 'listpix')

             ! cumul of unsorted tallies
             allocate(cumulu(0:nholes), stat=status)
             call assert_alloc(status, code, 'cumulu')
             cumulu(0) = 0
             do h=1,nholes
                cumulu(h) = cumulu(h-1) + tally(h-1)
             enddo
             ! cumul of sorted tallies
             allocate(cumuls(0:nholes), stat=status)
             call assert_alloc(status, code, 'cumuls')
             cumuls(0) = 0
             do h=1,nholes
                cumuls(h) = cumuls(h-1) + tally(indx(h-1))
             enddo

             do h=0,nholes-1
                ! define list boundaries
                k2 = cumuls(h)+nholes+1
                listpix(h) = k2

                ! fill list between boundaries
                np = cumuls(h+1) - cumuls(h)
                k = indx(h)
                listpix(k2:k2+np-1) = list(cumulu(k):cumulu(k+1)-1)
             enddo
             ! make sure upper boundary is correct
             listpix(nholes) = nph+nholes+1

             deallocate(cumuls, cumulu)
          endif
          deallocate(indx)
       endif
    endif

    deallocate(tally)
    deallocate(list)

!     ! revert Mask to its initial values
!     where( mask < 0) mask = 0

    return

  end subroutine size_holes_nest

end module mask_tools

