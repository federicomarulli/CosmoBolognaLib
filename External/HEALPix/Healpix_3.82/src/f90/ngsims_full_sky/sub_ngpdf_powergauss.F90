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
module sub_ngpdf_powergauss

  private
  public :: powergauss_driver
contains

!-----*-----------------------------------------------------------------

  SUBROUTINE powergauss_driver(handle, ns, sigma0, x, mu, bins)

    ! Draws a white noise map realisation from the 2nth power of a Gaussian
    ! variable with standard deviation sigma0. Its centred moments are mu
    ! (returns mu_1 ... mu_3). These are not the same as the cumulants for
    ! mu_4 and above.
    
    Use healpix_types
    Use paramfile_io, ONLY : paramfile_handle, parse_int, concatnl, parse_lgt
    use rngmod,    only: rand_init, planck_rng
    Use sky_sub
    
    implicit none
    integer(I4B), parameter :: KMAP = SP
    type(paramfile_handle), Intent(InOut)    :: handle
    integer(I4B), Intent(In)           :: ns !,npixmax 
    Real(DP),   Intent(Out)                  :: sigma0
    Real(KMAP), Intent(Out), dimension(0:ns-1) :: x
    real(DP),   Intent(Out), dimension(1:3)   :: mu
    Integer(I4B), Intent(In), Optional :: bins

    integer(I4B)                      ::  i,iseed,npower
    real(DP)                          ::  skewrel
    real(DP)                          ::  term2,term4,term6
    real(KMAP)                        ::  xmin,xmax
    character(len=filenamelen)        :: description
    character(len=20)                 :: chline
    logical(LGT)                      :: plot, do_bin = .false.
#ifdef PGPLOT
    integer, parameter                ::  npts=200  
    real(SP), dimension(1:npts)       ::  xline,yline
    real(DP)                          ::  xstep,xval
    integer :: pgdev, pgopen
#endif
    character(len=*), parameter :: code = 'shodev_driver'
    type(planck_rng) :: rng_handle

    If (present(bins)) do_bin = .True.

    !      initialize variables
!    chline = 'idum'
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

!    idum=-iabs(idum)
    call rand_init(rng_handle, iseed)    

    sigma0=1.d0

    chline = 'npower'
    If (do_bin) call add_subscript(chline, bins)
    description = concatnl(&
         & "", &
         & " Input n (variate is then x^2n):")
    npower = parse_int(handle, chline, default=1, vmin =1, vmax = 4, descr=description)
    ! NB using a value greater than 4 seems to crash due to numerical overflow

    ! Calculate centred moments of distribution

    ! term2=(2n-1)!!, term4=(4n-1)!!, term6=(6n-1)!!

    term2=1.0_dp
    do i=1,2*npower-1,2
       term2=term2*real(i, DP)
    end do
    !      Write (*,*) 'Term 2 is ', term2

    term4=1.0_dp
    do i=1,4*npower-1,2
       term4=term4*real(i, DP)
    end do
    !      Write (*,*) 'Term 4 is ', term4

    term6=1.0_dp
    do i=1,6*npower-1,2
       term6=term6*real(i, DP)
    end do
    !      Write (*,*) 'Term 6 is ', term6

    mu(1)=0.0_dp
    mu(2)=(term4 - term2**2)*sigma0**(4*npower)
    mu(3)=(term6 - 3.0_dp*term2*term4 + 2.0_dp*term2**3)*sigma0**(6*npower)

    skewrel=mu(3)/mu(2)**(3.0_dp/2.0_dp)

    write(*,*) 'mu1=', mu(1)
    write(*,*) 'mu2=', mu(2)
    write(*,*) 'mu3=', mu(3)
    write(*,*) 'Relative skewness=', skewrel

    ! Generate the map of non-Gaussian white noise

    xmin=0_KMAP
    xmax=0_KMAP
    do i = 0, ns-1
!       x(i)=powergauss(idum,npower,sigma0,term2)
       x(i) = powergauss(rng_handle, npower, sigma0, term2)
       if (x(i).lt.xmin) xmin=x(i)
       if (x(i).gt.xmax) xmax=x(i)
    end do
    xmax=min(xmax, real(sqrt(mu(2)), kind=KMAP)  )

#ifdef PGPLOT
    ! Find out if plotting is required
    description = concatnl( "", &
         & " Plot histogram of pixel values? (True or False)")
    plot = parse_lgt(handle, 'plot', default=.false., descr=description)
    If (plot) then 
       !     Write (*,*) 'Plot histogram of pixel values...'
!        call pgbegin(0,'?',1,1)
111    pgdev = pgopen('?')
       if (pgdev < 0) goto 111
       call pghist(ns,x,xmin,xmax,200,0)
       call pglab('Value','Number of pixels','Non-Gaussian PowGauss White Noise Map')

       xstep=(xmax-xmin)/real(npts-1, DP)
       do i=1,npts
          xval=xmin+real(i-1, DP)*xstep+1.e-5
          xline(i)=xval
          yline(i)=pdf(xval,sigma0,npower,term2)
       end do
       yline=yline*real(ns)*(xmax-xmin)/200.0
       call pgsci(12)
       call pgline(npts,xline,yline)

       call pgend

    End If !plot
#endif

  END SUBROUTINE powergauss_driver

  !-----*-----------------------------------------------------------------

  ! Returns the non-Gaussian pdf, for plotting only

  REAL FUNCTION pdf(x,sigma0,npower,term2)

    use healpix_types
    implicit none
    integer(i4b)            :: npower
    real(DP)               :: x,sigma0,term2

    real(DP)               :: xp,xg

    ! Add in mean again
    xp = x + term2*sigma0**(2*npower)

    ! Get Gaussian variable
    xg=xp**(0.5d0/real(npower, DP))

    ! Gaussian PDF
    pdf = exp(-0.5d0*(xg/sigma0)**2) / ( sqrt(2.d0*pi)*sigma0 )

    ! Apply Jacobian
    pdf = pdf /(real(npower, DP)*xg**(2.d0*npower-1.d0)) 

    return
  END FUNCTION pdf


  !-----*-----------------------------------------------------------------

  ! Draws a single sample from a Gaussian with standard deviation sigma0,
  ! raises to the power 2*npower and returns this as the function value

!  REAL FUNCTION powergauss(idum,npower,sigma0,term2)
  FUNCTION powergauss(rng_handle, npower, sigma0, term2) result (pgauss)

    use healpix_types
    use rngmod, only: planck_rng, rand_gauss

    implicit none
    real(DP) :: pgauss
    type(planck_rng), intent(inout) :: rng_handle
!    integer(i4b), intent(in)  :: idum
    integer(i4b), intent(in)  :: npower
    real(DP), intent(in)      :: sigma0,term2

    real(DP)               :: normal, cvar

    ! Get a Gaussian variate with unit variance

    normal = rand_gauss(rng_handle)

    ! Raise to 2*npower and centre
    cvar = normal**(2*npower)-term2

    ! Normalise
    pgauss = cvar*sigma0**(2*npower)

    return
  END FUNCTION powergauss


  !-----*-----------------------------------------------------------------

end module sub_ngpdf_powergauss























