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
module sub_ngpdf_sho


  use healpix_types
  integer(i4b), parameter :: namax = 20
  real(DP), dimension(0:namax) :: factorial
!  logical(LGT) :: precomp_done = .false.
  logical(LGT), save :: precomp_done ! 2013-05-07: edited for G95


  private
  public :: shodev_driver
  ! subroutine shodev_driver
  ! subroutine compute_factorial
  ! function shopdf
  ! subroutine rand_sho

contains


  !-----*-----------------------------------------------------------------

  SUBROUTINE shodev_driver(handle, ns, sigma0, x, nu, bins)

    Use healpix_types
    Use paramfile_io, ONLY : paramfile_handle, parse_int, parse_double, concatnl, parse_lgt
    use rngmod,    only: rand_init, planck_rng, rand_uni
    use misc_utils, only : wall_clock_time, assert_alloc
    Use sky_sub
    use statistics, only: compute_statistics, print_statistics, tstats

    implicit none

    !       type(paramfile_handle), Intent(In) :: handle
    integer(I4B), parameter :: KMAP = SP
    type(paramfile_handle), Intent(InOut)       :: handle
    Integer(I4B), Intent(In) :: ns
    Real(DP),   Intent(Out) :: sigma0
    Real(KMAP), Intent(Out), dimension(0:ns-1) :: x
    Real(DP),   Intent(Out) :: nu(1:3)
    Integer,    Intent(In), Optional :: bins

    integer                           ::  namax,nsmax 
    integer                           ::  i,j,iseed,na,status
    real(KMAP)                            ::  xmin,xmax
    real(DP)                              ::  skewrel, alpha0, mean
    real(DP), dimension(:), allocatable   ::  alpha
    character(len=filenamelen)          :: description
    character(len=20)                  :: chline
    logical :: plot, do_bin = .false.
    character(len = 20) :: form1 = '(a, i1)', form2 = '(a, i1, a)'
#ifdef PGPLOT
    Integer, parameter                ::  npts = 1000 
    real(SP), dimension(1:npts)       ::  xline,yline
    real(DP)                          ::  xval,xstep
    integer :: pgopen, pgdev
#endif
    character(len=*), parameter :: code = 'shodev_driver'
    type(planck_rng) :: rng_handle
    real(DP), parameter :: one = 1.0_dp, two = 2.0_dp, three = 3.0_dp
    real(DP), parameter :: four = 4.0_dp, six = 6.0_dp
    real(DP), allocatable, dimension(:) :: shovar
    real :: time0, time1, time2
    integer :: loc(1)
    type(tstats) :: my_stats


    !      initialize variables

    If (present(bins)) do_bin = .True.
    precomp_done = .false.

    chline = 'iseed'
    If (do_bin) call add_subscript(chline, bins)
    description = concatnl( "", &
         & " Input seed integer iseed:")
    iseed = parse_int(handle, chline, default=1, descr=description)

    !     --- If necessary use current system time to set seed  ---
    IF (iseed .eq. 0) THEN
       CALL SYSTEM_CLOCK(COUNT = iseed)
       !-- Seed should initially be set to a large, odd integer --
       IF (MOD(iseed,2) .EQ. 0) iseed = iseed + 1
       PRINT *,"      "//code//"> Generating random number seed ", iseed
    END IF

    call rand_init(rng_handle, iseed)    

    chline = 'sigma0'
    If (do_bin) call add_subscript(chline, bins)
    description = concatnl( "", &
         & " Input value of sigma0:")
!    sigma0 = parse_real(handle, chline, default=1.0, descr=description)
    sigma0 = parse_double(handle, chline, default=1.d0, descr=description)

    chline = 'na'
    If (do_bin) call add_subscript(chline, bins)
    description = concatnl( "", &
         & " Input order of highest non-zero alpha (na):")
    na = parse_int(handle, chline, default=3, vmin = 0, vmax = 20, descr=description)
    ! Note maximum value of na is 20, but currently moments can only be calculated analytically
    ! for na no greater than 3


    !---allocate space for arrays---
    allocate(shovar(0:ns-1), stat=status)
    call assert_alloc(status, code, ' cannot allocate space for alpha array')

    allocate(alpha(1:MAX(na,3)), stat=status)
    call assert_alloc(status, code, ' cannot allocate space for alpha array')
    alpha = 0.0_dp
    !---end of allocation for arrays---

    Do i = 1, na
       If (i .EQ. 10) Then
          form1 = "(a, i2)"
          form2 = "(a, i2, a)"
       End If
       Write (chline, form1) "alpha_", i
       Write (description, form2) ' Input alpha(',i,')'
       If (do_bin) call add_subscript(chline, bins)
       description = concatnl( "", description)
       alpha(i) = parse_double(handle, chline, default = 0.0_dp, vmin = 0.0_dp, &
            &vmax = sqrt(1.0_dp-SUM(alpha**2)), descr = description)
    End Do

    ! 
    Write (*,*) 'Drawing pixel values from non-Gaussian pdf...'
    if (.not. precomp_done) call compute_factorial()
    call rand_sho(rng_handle, ns, shovar, sigma0, alpha, na)
    call compute_statistics(shovar, my_stats)
    call print_statistics(my_stats)
    xmin = minval(shovar)
    xmax = maxval(shovar)

    mean = 0.0
    if (na.le.3) then
       !      compute moments of this distribution
       alpha0 = sqrt(1.0 - SUM(alpha**2))
       write (*,*) 'alpha0 is ',alpha0
       nu(1)=sqrt(TWO*sigma0)*(TWO*alpha(1)*alpha(2) + sqrt(TWO)*alpha0*alpha(1) &
            & + sqrt(SIX)*alpha(2)*alpha(3))
       nu(2)=sigma0**2 *(ONE + TWO*alpha(1)**2 + FOUR*alpha(2)**2 + SIX*alpha(3)**2 &
            & + TWO*sqrt(TWO)*alpha0*alpha(2) + TWO*sqrt(SIX)*alpha(1)*alpha(3))
       nu(3)=THREE*alpha0*alpha(1)/sqrt(TWO) + SIX*alpha(1)*alpha(2) + &
            & 9.0_dp*sqrt(THREE)*alpha(2)*alpha(3)/sqrt(TWO) + sqrt(THREE)*alpha0*alpha(3)
       nu(3)=(TWO*sigma0**2)**(THREE/TWO)*nu(3)

       skewrel=nu(3)/nu(2)**(THREE/TWO)

       Write(*,*) 'Raw (theoretical) moments of distribution:'
       write(*,*) 'nu1=', nu(1)
       write(*,*) 'nu2=', nu(2)
       write(*,*) 'nu3=', nu(3)
       write(*,*) 'Relative skewness=', skewrel

       If (nu(1) .NE. 0.0) Then
          Write(*,*) 'Re-centering distribution on zero'
          mean = nu(1)
          shovar = shovar - nu(1)
          xmin = xmin - nu(1)
          xmax = xmax - nu(1)
          nu(3) = nu(3) - THREE*nu(1)*nu(2) + TWO*nu(1)**3
          nu(2) = nu(2) - nu(1)**2
          nu(1) = 0.0_dp
          skewrel=nu(3)/nu(2)**(THREE/TWO)
          Write (*,*) 'Centered (theoretical) moments are:'
          Write (*,*) 'nu1=', nu(1)
          write (*,*) 'nu2=', nu(2)
          write (*,*) 'nu3=', nu(3)
          write (*,*) 'Relative skewness=', skewrel
       End If

    Else
       nu = 0.0_dp
       nu(1) = -1.0_dp
       Write (*,*) 'Cannot calculate moments when highest non-zero alpha > 3'
    End If

    x = shovar

    ! Find out if plotting is required
#ifdef PGPLOT
    description = concatnl( "", &
         & " Plot histogram of pixel values? (True or False)")
    plot = parse_lgt(handle, 'plot', default=.false., descr=description)

    If (plot) then
!        call pgbegin(0,'?',1,1)
111    pgdev = pgopen('?')
       if (pgdev < 0) goto 111
       call pghist(ns,real(x,SP),real(xmin,SP),real(xmax,SP),200,0)
       call pglab('Pixel value','Number of pixels','Non-Gaussian SHO White Noise Map')

       xstep=(xmax-xmin)/real(npts-1, DP)
       Do i=1,npts        
          xval=xmin+real(i -1, DP)*xstep
          xline(i)=xval
          yline(i)=shopdf(xval+mean, sigma0, alpha, na)
       End Do
       yline = yline*real(ns)*(xmax-xmin)/200.
       call pgsci(12)
       call pgline(npts,xline,yline)
       call pgend
    End If !plot
#endif


  END SUBROUTINE shodev_driver

  !==================================================================================
  subroutine compute_factorial()
    !-------------------------------------------------------------------------------
    integer(i4b) :: i


    if (precomp_done) return

    factorial(0) = 1.0_dp
    factorial(1) = 1.0_dp
    
    do i = 2, namax
       factorial(i) = factorial(i-1) * real(i, DP)
    enddo
    precomp_done = .true.

    return
  end subroutine compute_factorial

  !-----*-----------------------------------------------------------------

  FUNCTION shopdf(x, sigma0, alpha, na) result(spdf)
    !------------------------------------------------------------
    ! calculates the general simple harmonic oscillator pdf 
    !------------------------------------------------------------
    use healpix_types
    use num_rec, only : othpl
    implicit none
    real(DP), intent(in)                 :: x
    integer(i4b), intent(in)             :: na
    real(DP), intent(in)                 :: sigma0
    real(DP), intent(in)                 :: alpha(Max(na,3))
    real(DP)                             :: spdf

    integer, parameter   :: namax=20
    integer(i4b)         :: i,m,n
    real(DP), dimension(0:namax)     :: pars, hermxp, dhermxp
    real(DP)                         :: xp,sumsq,fact
    !------------------------------------------------------------

    xp = x/(sqrt(2.d0)*sigma0)

    pars(1:na) = alpha(1:na)
    sumsq=0d0
    do i=1,na
       sumsq=sumsq + pars(i)*pars(i)
    end do
    pars(0)=sqrt(1d0-sumsq)

    call othpl(4,na,xp,hermxp,dhermxp)

    if (.not. precomp_done) call compute_factorial

    spdf=0.0_dp
    do m=0,na
       do n=0,na
          fact=  sqrt( (2d0**dble(m+n)) * factorial(m) * factorial(n) )
          spdf=spdf+ pars(m)*pars(n)*hermxp(m)*hermxp(n) / fact
       end do
    end do
    spdf=spdf* exp(-xp*xp)/sqrt(2d0*pi)/sigma0

    return
  END FUNCTION shopdf

! =====================================================================
  subroutine rand_sho(rng_handle, npix, xvar, sigma0, alpha, na)
    !------------------------------------------------------------
    ! draws npix numbers from a SHO variable
    !------------------------------------------------------------
    use healpix_types
    use misc_utils, only: fatal_error
    use rngmod, only: rand_uni, planck_rng
    implicit none

    type(planck_rng), intent(inout)      :: rng_handle
    integer(i4b), intent(in)             :: na, npix
    real(DP), intent(out)                :: xvar(0: npix-1)
    real(DP), intent(in)                 :: sigma0
    real(DP), intent(in)                 :: alpha(Max(na,3))

    character(len=*), parameter :: code = 'rand_sho'
    real(DP), parameter :: threshold = 1.d-12

    real(DP) :: x, xmin, xmax, xstep, xlin, y, yold, yuni
    real(DP), dimension(:), allocatable :: integral
    integer(I4B) :: nx, i, p, ilow, ihi, imid
    !----------------------------------------------------------

    ! find out x range
    xmax =  5.0_dp * sigma0
    xmin = -xmax
    do
       if (shopdf(xmax, sigma0, alpha, na) < threshold) exit
       xmax = xmax + sigma0
    enddo
    do
       if (shopdf(xmin, sigma0, alpha, na) < threshold) exit
       xmin = xmin - sigma0
    enddo
    nx = (xmax - xmin)/sigma0 * 2000

    ! build cumulative table, using trapeze integral
    allocate(integral(0:nx-1))
    xstep = (xmax-xmin)/(nx-1.d0)
    integral = 0_dp
    yold = 0.0_dp
    do i=0, nx - 1
       x = xmin +  i * xstep
       y = shopdf(x, sigma0, alpha, na)
       if (i > 0) integral(i) = integral(i-1) + (y+yold)*0.5_dp
       yold = y
    enddo
    integral = integral * xstep ! normalize
    if (abs(integral(nx-1)-1.d0) > 5.d0*threshold) then
       call fatal_error(code//': problem with computation of integral')
       print*,'Error on integral: ',integral(nx-1)-1.d0
    endif
    integral(nx-1) = 1.d0
    
    ! solve y = F(x)
    do p = 0, npix-1
       ! uniform variate in [0,1]
       yuni = rand_uni(rng_handle)
       ! locate point in tabulated integral by dichotomy
       ilow = 0
       ihi  = nx-1
       do 
          imid = (ihi + ilow)/2
          if ((integral(imid) <= yuni) .eqv. (yuni <= integral(ihi))) then
             ilow = imid
          else
             ihi = imid
          endif
          if ((ihi - ilow) == 1) exit
       enddo
       ! linear interpolation
       xlin = (yuni - integral(ilow))/(integral(ihi)-integral(ilow)) * xstep
       ! return non-gaussian variate
       xvar(p) =  xmin + ilow * xstep + xlin
    enddo
    deallocate(integral)

    return
  end subroutine rand_sho
  !-----*-----------------------------------------------------------------

end module sub_ngpdf_sho






















