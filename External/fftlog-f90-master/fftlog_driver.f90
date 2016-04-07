! Driver for FFTLog-f90.
!
! FFTLog-f90 is the Fortran 90 version of FFTLog, a Fortran 77 code by
! Andrew Hamilton to compute the Fourier-Bessel, or Hankel, transform of
! of a series of logarithmically spaced points.
!
! The main reference for the algorithm is Andrew Hamilton's webpage,
! http://casa.colorado.edu/~ajsh/FFTLog. He first introduced FFTLog in a
! paper on the nonlinearities in the cosmological structures
! (http://xxx.lanl.gov/abs/astro-ph/9905191).
!
! The Fourier-Bessel, or Hankel, transform F(r) of a function f(k) is defined
! as in eq. 159 of Hamilton's paper:
!
!           /
!  F(r) =  | dk * k * J_mu(k*r) * f(k)
!          /
!
! where J_mu(k*r) is the Bessel function of order mu and argument k*r. 
!
! IMPORTANT: Currently, FFTLog-f90 only supports the sine transform (mu=0.5)
! and returns the following integral:
!
!              1       /            sin(k*r)
!  F(r) = ----------  | dk * k^2 * ---------- * f(k)
!          (2*pi)^2   /                k*r
!
! If f(k)=P(k) is the power spectrum of a homogeneous 3D random field, then the
! integral will yield the two-point correlation function xi(r), as in eq. 3.104
! of http://arxiv.org/abs/1405.2280.
!
! To generalise FFTLog-f90 to arbitrary Bessel order, one needs to use
! the FHTQ subroutine rather than the FFTL one.
!
! ARGUMENTS
!
! This is the list of command-line arguments taken by FFTLog-f90:
!
! 1. The path of the input file, containing the k and f(k) columns.
!
! 2. The path of the output file, containing the r and F(r) columns. The file
!    will have as many rows as the input file, unless a third argument is given.
!    By default, the range in r will go from 1/k_max to 1/k_min and will contain
!    as many points as the rows in the input file.
!
! 3. [OPTIONAL] N_OUTFILE, the number of r values to include in the output file.
!    Set to zero to use the same number of points in the input file (default behaviour).
!    If positive, then f(k) will be interpolated in N_OUTFILE points and fed to FFTLog.
!
! 4. [OPTIONAL] The order mu of the Bessel function J_mu. By default we take mu=0.5,
!    which corresponds to a regular Fourier transform.
!
! 5. [OPTIONAL] The bias q of the Hankel transformation, as defined in eq. 156.
!    By default we take it to vanish (b=0).
!
! 6. [OPTIONAL] Inferior limit in k for the integration. By default it is set 
!    to the first value in the k column.
!
! 7. [OPTIONAL] Superior limit in k for the integration. By default it is set 
!    to the first value in the k column.
!
! Created by Guido Walter Pettinari on 26/07/2010
! Last modified on 09/09/2015 by GWP


PROGRAM fftlog_f90

	IMPLICIT NONE

	REAL ( KIND = 8 ), PARAMETER :: PI=3.141592653589793238462643383279502884197d0, CONSTANT= 1/(2*PI*PI)*sqrt(PI/2)
	REAL ( KIND = 8 ) :: xx_inf_limit = -1d300, xx_sup_limit = 1d300, dln_xx, dlog_xx, log_xmedian,&
                       &log_xxmedian, log_xxmax, log_xxmin, rk, kr, mu, CENTRAL_INDEX, q, x, temp1, temp2
	INTEGER :: N_ARGUMENTS, DIRECTION, KR_OPTION, UNIT, error, N_INFILE=0, N_USED=0, FIRST_USED_INDEX = 0,&
             &LAST_USED_INDEX=0, N_OUTFILE=0, i=0
	LOGICAL :: INTERPOLATION = .FALSE., FFT_OK
	CHARACTER(512) :: input_filename, output_filename
	CHARACTER(128) :: buffer
  !	The arrays xx and yy will contain the x and y values from the input file, respectively.
  ! They will be overwritten by the FFTL subroutine with the Hankel transform.
	REAL ( KIND = 8 ), DIMENSION(:), ALLOCATABLE :: xx, yy
  ! Parameters for the spline interpolation
	REAL ( KIND = 8 ), DIMENSION(:), ALLOCATABLE :: yy_second_derivative, yy_interp, wsave
  !	Complain if the input file has more than N_MAX lines, because it would require
  ! about 2 GB or RAM. Increase N_MAX if you have more memory than that.
	INTEGER, PARAMETER :: N_MAX=67108864
	
  
  
  ! =======================================================================================
  ! =                                   Parse arguments                                   =
  ! =======================================================================================
  
	N_ARGUMENTS = COMMAND_ARGUMENT_COUNT()
	
	IF ( N_ARGUMENTS.GT.7 .OR. N_ARGUMENTS.LT.2 ) STOP &
   & "FFTLog-f95 can have between 2 and 7 arguments. For more details, please refer to the&
   & documentation in fftlog_driver.f90 (https://github.com/coccoinomane/fftlog-f90)."


	CALL GET_COMMAND_ARGUMENT( 1, input_filename )
	CALL GET_COMMAND_ARGUMENT( 2, output_filename )

	IF ( N_ARGUMENTS .GE. 3 ) THEN
		CALL GET_COMMAND_ARGUMENT( 3, buffer )
		READ ( buffer, * ) N_OUTFILE
	END IF

	! order of Bessel function ( mu = 0.5 for sine transform )
  mu = 0.5d0
	IF ( N_ARGUMENTS .GE. 4 ) THEN
		CALL GET_COMMAND_ARGUMENT( 4, buffer )
		READ ( buffer, * ) mu
	END IF
  IF (mu .NE. 0.5d0) STOP &
   & "FFTLog-f90 so far only supports sine (mu=0.5) and cosine (mu=0.5) transforms. Make&
   & sure that the fourth argument is either 0.5 or -0.5, respectively."

	! bias exponent: q = 0 is unbiased ( good for power spectrum <-> correlation function )
  q = 0.d0
	IF ( N_ARGUMENTS .GE. 5 ) THEN
		CALL GET_COMMAND_ARGUMENT( 5, buffer )
		READ ( buffer, * ) q
	END IF

	IF ( N_ARGUMENTS == 6 ) THEN
		CALL GET_COMMAND_ARGUMENT( 6, buffer )
		READ ( buffer, * ) xx_sup_limit
	END IF
		
	IF ( N_ARGUMENTS == 7 ) THEN
		CALL GET_COMMAND_ARGUMENT( 6, buffer )
		READ ( buffer, * ) xx_inf_limit
		CALL GET_COMMAND_ARGUMENT( 7, buffer )
		READ ( buffer, * ) xx_sup_limit
	END IF
	
  
  
  ! =======================================================================================
  ! =                                   Read input file                                   =
  ! =======================================================================================  
  
  ! Select the entries of input_filename to process
	OPEN ( UNIT=50, FILE=input_filename, STATUS='OLD', IOSTAT=error )
	IF ( error /= 0 ) STOP "Input file could not be opened, exiting..."
	DO
		READ (UNIT=50, FMT=*, IOSTAT=error) temp1
		IF ( error < 0 ) EXIT
		N_INFILE = N_INFILE + 1		
		IF ( temp1 < xx_inf_limit .OR. temp1 > xx_sup_limit ) CYCLE
		N_USED = N_USED + 1
		LAST_USED_INDEX = N_INFILE
	END DO
	CLOSE(UNIT = 50)
	FIRST_USED_INDEX = LAST_USED_INDEX - N_USED + 1
	
  !	Control block	
	IF ( N_USED == 0 ) STOP "Either the input file is empty or incompatible with the specified integration limits."	
	IF ( N_OUTFILE .LE. 0 ) THEN	
		N_OUTFILE = N_USED
		INTERPOLATION = .FALSE.
!		PRINT *, "No interpolation of input data will be performed."
	ELSE
		INTERPOLATION = .TRUE.
	END IF		
    IF ( (N_USED > N_MAX) .OR. (N_OUTFILE > N_MAX) ) &
    &STOP "The number you specified or the number of elements of the input file are greater that N_MAX"

  ! It is now safe to allocate memory to the arrays
	ALLOCATE( xx(N_USED), yy(N_USED), yy_second_derivative(N_USED) )
  ALLOCATE( yy_interp(N_OUTFILE), wsave ( 2*N_OUTFILE + 3*(N_OUTFILE/2) + 19 ) )

  ! Read the first two columns of the input file
	OPEN ( UNIT=50, FILE=input_filename, STATUS='OLD', IOSTAT=error )
	IF ( error /= 0 ) STOP "Input file could not be opened, exiting..."
  !	Dump the rows of the file that are smaller than xx_inf_limit
	DO i = 1, FIRST_USED_INDEX - 1
		READ ( 50, FMT=*, IOSTAT=error )
	END DO
	DO i = 1, N_USED
		READ ( 50, FMT=*, IOSTAT=error ) temp1, temp2
		IF ( error < 0 ) EXIT
		xx(i) = temp1
		yy(i) = temp2
	END DO
	CLOSE( UNIT = 50 )
	
	!PRINT *, "Inferior integration limit = ", xx(1)
	!PRINT *, "corresponding in input file to row = ", FIRST_USED_INDEX
	!PRINT *, "Superior integration limit = ", xx(N_USED)
	!PRINT *, "corresponding in input file to row = ", LAST_USED_INDEX	
	!PRINT *, "Order of Bessel function mu = ", mu
	!PRINT *, "Number of elements in input file = ", N_INFILE
  !PRINT *, "Number of elements used = ", N_USED
	!PRINT *, "Number of elements in output file = ", N_OUTFILE



  ! =======================================================================================
  ! =                                Prepare integrand f(k)                               =
  ! =======================================================================================  

  ! Take the logarithm of the x limits
	log_xxmin = LOG10( xx(1) )
	log_xxmax = LOG10( xx(N_USED) )
  
  ! Logarithmic step in x
	dlog_xx = (log_xxmax-log_xxmin) / (N_OUTFILE-1)

	! central index (1/2 integral if N_OUTFILE is even)
	CENTRAL_INDEX = dble( N_OUTFILE+1 ) / 2.d0
	
	! logarithmical spacing between points (needed by fhti). Was dlnr
	dln_xx = dlog_xx * LOG(10.d0)

	! central point of periodic interval at log10 (xxmedian). Was logrc
	log_xxmedian = ( log_xxmin + log_xxmax) / 2.d0

	! sensible approximate choice of k_c r_c
	kr = 1.d0

	! tell fhti to change kr to low-ringing value
	KR_OPTION = 1

	! forward transform
	DIRECTION = 1

	! central point in k-space. Was logkc
	log_xmedian = LOG10 (kr) - log_xxmedian

	! rk = r_c/k_c
  rk = 10.d0 ** ( log_xxmedian - log_xmedian )

  ! Interpolate the input 
	IF( INTERPOLATION ) THEN
		! Compute the second derivative of the input array yy. The prototype is
    ! spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp )
		CALL SPLINE_CUBIC_SET ( N_USED, xx, yy, 0, 0.0d0, 0, 0.0d0, yy_second_derivative )
		DO i = 1, N_OUTFILE ! the x's are the input domain
			x = 10.d0 ** ( log_xxmedian + ( i - CENTRAL_INDEX )*dlog_xx )
      ! Interpolate the input array yy in x. The prototype is
			! spline_cubic_val ( n, t, y, ypp, tval, yval, ypval, yppval )
			CALL SPLINE_CUBIC_VAL ( N_USED, xx, yy, yy_second_derivative, x, yy_interp(i), temp1, temp2 )
			yy_interp(i) = yy_interp(i)*x
		END DO
		DEALLOCATE ( yy_second_derivative )
	ELSE
		DO i = 1, N_USED
			yy(i) = xx(i) * yy(i)
		END DO
	END IF
  
  
  
  ! =======================================================================================
  ! =                                    Call FFTLog                                      =
  ! =======================================================================================

  ! Initialize FFTLog transform. Note that fhti resets kr
	CALL FHTI ( N_OUTFILE, mu, q, dln_xx, kr, KR_OPTION, wsave, FFT_OK )
	IF( .NOT. FFT_OK ) STOP "FHTI not ok!"

  ! Call FFTLog
	IF( INTERPOLATION ) THEN 
		CALL FFTL ( N_OUTFILE, yy_interp, rk, DIRECTION, wsave )
	ELSE
		CALL FFTL ( N_OUTFILE, yy, rk, DIRECTION, wsave )
	END IF
	


  ! =======================================================================================
  ! =                                Write result to file                                 =
  ! =======================================================================================

	OPEN( UNIT=40, FILE = output_filename, STATUS = "REPLACE", IOSTAT = error )
	IF ( error /= 0 ) STOP "Output file could not be opened, exiting..."
	IF ( INTERPOLATION ) THEN 
		DO i = 1, N_OUTFILE ! now the x's are the output domain
			x = 10.d0 ** ( log_xmedian + ( i - CENTRAL_INDEX ) * dlog_xx )
			WRITE ( 40, '(3g24.16)' ) x, CONSTANT * yy_interp(i)/x
		END DO
		DEALLOCATE ( yy_interp )
	ELSE
		DO i = 1, N_OUTFILE ! now the x's are the output domain
			x = 10.d0 ** ( log_xmedian + ( i - CENTRAL_INDEX )*dlog_xx )
			WRITE( 40, '(3g24.16)' ) x, CONSTANT * yy(i)/x
		END DO
	END IF
	DEALLOCATE(xx,yy,wsave)
	
END PROGRAM fftlog_f90


