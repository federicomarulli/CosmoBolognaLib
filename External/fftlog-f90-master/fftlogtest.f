c-----------------------------------------------------------------------
c * This is fftlogtest.f
c *
c * This is a simple test program to illustrate how FFTLog works.
c * The test transform is
c *    inf
c *    /   mu+1                               mu+1
c *    |  r     exp(-r^2/2) J_mu(k r) k dr = k     exp(-k^2/2)
c *    /
c *   0
c *
c * Disclaimer:
c * FFTLog does NOT claim to provide the most accurate possible
c * solution of the continuous transform (which is the stated aim
c * of some other codes).  Rather, FFTLog claims to solve the exact
c * discrete transform of a logarithmically-spaced periodic sequence.
c * If the periodic interval is wide enough, the resolution high
c * enough, and the function well enough behaved outside the periodic
c * interval, then FFTLog may yield a satisfactory approximation
c * to the continuous transform.
c *
c * Observe:
c * (1) How the result improves as the periodic interval is enlarged.
c *     With the normal FFT, one is not used to ranges orders of
c *     magnitude wide, but this is how FFTLog prefers it.
c * (2) How the result improves as the resolution is increased.
c *     Because the function is rather smooth, modest resolution
c *     actually works quite well here.
c * (3) That the central part of the transform is more reliable
c *     than the outer parts.  Experience suggests that a good general
c *     strategy is to double the periodic interval over which the
c *     input function is defined, and then to discard the outer
c *     half of the transform.
c * (4) That the best bias exponent seems to be q = 0.
c * (5) That for the critical index mu = -1, the result seems to be
c *     offset by a constant from the `correct' answer.
c * (6) That the result grows progressively worse as mu decreases
c *     below -1.
c *
c * The analytic integral above fails for mu <= -1, but FFTLog
c * still returns answers.  Namely, FFTLog returns the analytic
c * continuation of the discrete transform.  Because of ambiguity
c * in the path of integration around poles, this analytic continuation
c * is liable to differ, for mu <= -1, by a constant from the `correct'
c * continuation given by the above equation.
c *
c * FFTLog begins to have serious difficulties with aliasing as
c * mu decreases below -1, because then r^(mu+1) exp(-r^2/2) is
c * far from resembling a periodic function.
c * You might have thought that it would help to introduce a bias
c * exponent q = mu, or perhaps q = mu+1, or more, to make the
c * function a(r) = A(r) r^-q input to fhtq more nearly periodic.
c * In practice a nonzero q makes things worse.
c *
c * A symmetry argument lends support to the notion that the best
c * exponent here should be q = 0, as empirically appears to be true.
c * The symmetry argument is that the function r^(mu+1) exp(-r^2/2)
c * happens to be the same as its transform k^(mu+1) exp(-k^2/2).
c * If the best bias exponent were q in the forward transform, then
c * the best exponent would be -q that in the backward transform;
c * but the two transforms happen to be the same in this case,
c * suggesting q = -q, hence q = 0.
c *
c * This example illustrates that you cannot always tell just by
c * looking at a function what the best bias exponent q should be.
c * You also have to look at its transform.  The best exponent q is,
c * in a sense, the one that makes both the function and its transform
c * look most nearly periodic.
c *
c        parameters
      integer NMAX
      parameter (NMAX=4096)
c        externals
      integer lnblnk
c        local variables
      character*128 outfile
      integer dir,i,kropt,n,unit
      logical ok
      real*8 a(NMAX),dlnr,dlogr,k,kr,logkc,logrc,logrmax,logrmin,
     *  mu,nc,q,r,rk
      real*8 wsave(2*NMAX+3*(NMAX/2)+19)
c
      print *,'test integral_0^inf r^(mu+1) exp(-r^2/2) J_mu(k r) k dr =
     * k^(mu+1) exp(-k^2/2)'
c--------reasonable choices of parameters
c        range of periodic interval
      logrmin=-4.d0
      logrmax=4.d0
c        number of points
      n=64
c        order of Bessel function
      mu=0.d0
c        bias exponent: q = 0 is unbiased
      q=0.d0
c        sensible approximate choice of k_c r_c
      kr=1.d0
c        tell fhti to change kr to low-ringing value
      kropt=3
c        forward transform
      dir=1
c--------choose parameters
  100 print *,'enter period range log10(r_min), log10(r_max) [',
     *  logrmin,',',logrmax,']'
      read (*,*,end=300,err=300) logrmin,logrmax
c        central point log10(r_c) of periodic interval
      logrc=(logrmin+logrmax)/2.d0
      print *,'central point of periodic interval at log10(r_c) =',logrc
      print *,'enter number of points [',n,']'
      read (*,*,end=300,err=300) n
      if (n.gt.NMAX) then
        print *,n,' exceeds declared maximum',NMAX
        goto 100
      endif
c        central index (1/2 integral if n is even)
      nc=dble(n+1)/2.d0
c        log spacing of points
      dlogr=(logrmax-logrmin)/n
      dlnr=dlogr*log(10.d0)
c        order mu of Bessel function
      print *,'enter order mu of Bessel function [',mu,']'
      read (*,*,end=300,err=300) mu
c        bias exponent q
      print *,'enter bias exponent q [',q,']'
      read (*,*,end=300,err=300) q
c        kr = k_c r_c
      print *,'enter kr = k_c r_c [',kr,']'
      read (*,*,end=300,err=300) kr
c--------r^(mu+1) exp(-r^2/2)
      do i=1,n
        r=10.d0**(logrc+(i-nc)*dlogr)
        a(i)=r**(mu+1.d0)*exp(-r**2/2.d0)
      enddo
c--------initialize FFTLog transform - note fhti resets kr
      call fhti(n,mu,q,dlnr,kr,kropt,wsave,ok)
      if (.not.ok) goto 100
      logkc=log10(kr)-logrc
      print *,'central point in k-space at log10(k_c) =',logkc
c        rk = r_c/k_c
      rk=10.d0**(logrc-logkc)
c--------transform
c     call fftl(n,a,rk,dir,wsave)
      call fht(n,a,dir,wsave)
c--------print/write result
      print *,'enter name of file to save transform [CR=screen]'
      read (*,'(a128)') outfile
      if (outfile.eq.' ') then
        unit=6
      else
        unit=7
        open (unit,file=outfile)
        rewind (unit)
      endif
      write (unit,'(3a24)') 'k    ','a(k)    ','k^(mu+1) exp(-k^2/2)'
      do i=1,n
        k=10.d0**(logkc+(i-nc)*dlogr)
        write (unit,'(3g24.16)') k,a(i),k**(mu+1.d0)*exp(-k**2/2.d0)
      enddo
      call flush(unit)
      if (outfile.ne.' ') then
        print *,'header +',n,' lines written to file ',
     *    outfile(1:lnblnk(outfile))
      endif
c--------end of test loop
      goto 100
  300 continue
      end
c
