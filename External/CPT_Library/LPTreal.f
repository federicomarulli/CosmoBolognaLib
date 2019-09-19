c *********************************************************  c
c                                                            c
c      Lagrangian Perturbation Theory for                    c
c           Non-linver Evolution of Power Spectrum           c
c                                                            c
c             Time-stamp: <2010-09-18 15:41:57 ataruya>      c
c                                                            c
c     In addition to the standard PT power spectra, pk22     c
c     and pk13, this codes calculate the term, pk_corr,      c  
c     which can be used to compute power spectrum in         c
c     Lagrangial PT developed by Matsubara (arXiv:0711.2521) c
c                                                            c
c     Note--.                                                c
c                                                            c
c     The correction pk_corr can be used to compute P(k)     c
c     for Lagrangian PT as                                   c
c                                                            c
c     P(k) = exp(- D^2 * pk_corr )                           c
c                     * ( D^2*pk11 + D^4*(pk22 + pk13) )     c
c                                                            c
c     Notice that this code only provides the data for PT    c
c     corrections,  i.e., pk13, pk22, pk_corr. To compute    c
c     non-linear P(k) at arbitrary redshift, you should run  c
c     the common read file, read_pk.f                        c
c                                                            c
c ********************************************************** c
c
c     Note--. 
c
c     All the parameter 'ikmax' in this program must be 
c     the same. 
c
      program LPTreal
c
      implicit none
c
      integer  ik, ikmax, ik_max
      parameter(ikmax=3000)
      real*8   ak(ikmax), pk(ikmax), pk22(ikmax), pk13(ikmax)
      real*8   pk_corr(ikmax)
      real*8   sigma8, kmin, kmax
      common /normalized_pk/ ak, pk, ik_max
c     -------------------------------------------------
c
      write(6,*)
      write(6,*) '*** Lagrangian PT calculation for power spectrum ***'
      write(6,*)
c
      write(6,*) 'Input sigma_8:' 
      read(5,*) sigma8
c
c     /////// Load transfer function ///////
c
      call load_transfer_function
c
      write(6,*) ' loading (linear) matter power spectrum, done '
c
c     /////// Sigma8 normalization ///////
c
      call normalization_trapez(sigma8)
c
      write(6,*) ' normalization done '
c
c     /////// range of wavenumber used for 1-loop correction ///////
c     
c     ! default values 
      kmin = 5.d-4
      kmax = 10.d0
c     
      call truncation_k(kmin, kmax)
c
c     /////// 1-loop correction of P(k) ///////
c
      call calc_one_loop_pk(pk13, pk22) 
c
      write(6,*) ' 1-loop P(k) done '
c
c     /////// Corrections by T.Matsubara (arXiv:0711.2521) ///////
c
      call calc_correction(pk_corr)
c
      write(6,*) ' Resummed correction done '
c
c     /////// Save output data ///////
c
      open(10,file='LPT.dat',status='unknown')
c
      do ik =1, ik_max
         write(10,'(1p6e18.10)') 
     &        ak(ik), pk(ik), pk13(ik), pk22(ik), pk_corr(ik)
      enddo
c
      close(10)
c
      end
c
c ******************************************************* c
c
      subroutine load_transfer_function
c
c ******************************************************* c
c
c     input file is assumed to be the matter power spectrum data 
c     created by CAMB code, that is, the data sorting is assumed
c     to be 
c
c     k [h/Mpc],  P(k) [Mpc^3/h^3]
c
      implicit none
c
      integer ik, ikmax, ik_max
      parameter(ikmax=3000)
      real*8  ak(ikmax), pk(ikmax)
      common /normalized_pk/ ak, pk, ik_max
      character infile*50
c     --------------------------------------------------------
c
      write(6,*) ' Type file name of matter P(k) data'
      read(5,*) infile
      open(9, file=infile, status='unknown')
c
      do ik=1, ikmax
         read(9,*,END=10) ak(ik), pk(ik)
      enddo
c
 10   continue
      ik_max = ik - 1
c
      write(6,*) 'ak(1)=', ak(1)
      write(6,*) 'ak(ik_max)=', ak(ik_max)
      write(6,*) 'ik_max=', ik_max
c
      close(9)
c
      end
c
c ******************************************************* c
c
      subroutine normalization_trapez(sigma8)
c
c ******************************************************* c
c
      implicit none
      integer ik, ik_max, ikmax
      parameter(ikmax=3000)
      real*8  sigma8
      real*8  ak(ikmax), pk(ikmax)
      real*8  r_th, const, x
      real*8  W_TH, sigma_a, sigma_b, pi
      common /normalized_pk/ ak, pk, ik_max
      pi = 4.d0 * atan(1.d0)
      r_th = 8.d0
c     ---------------------------------------------------
c
      x = ak(1) * r_th
      if(x.lt.1.d-3) then
         W_TH = 1.d0 - x*x / 10.d0 + x**4 / 280.d0 
      else
         W_TH = 3.d0 * (sin(x) - x * cos(x))/x/x/x
      endif
      sigma_a = W_TH * W_TH * pk(1) * ak(1) * ak(1)
      sigma_a = sigma_a / (2.d0 * pi * pi)
c
      const = 0.d0 
      do ik=2, ik_max
         x = ak(ik) * r_th
         if(x.lt.1.d-3) then
            W_TH = 1.d0 - x*x / 10.d0 + x**4 / 280.d0 
         else
            W_TH = 3.d0 * (sin(x) - x * cos(x))/x/x/x
         endif
         sigma_b = W_TH * W_TH * pk(ik) * ak(ik) * ak(ik) 
         sigma_b = sigma_b / (2.d0 * pi * pi)
         const = const + 
     &        (sigma_a + sigma_b) * ( ak(ik) - ak(ik-1) )/ 2.d0
         sigma_a = sigma_b
      enddo
c
      do ik=1, ik_max
         pk(ik) = sigma8 * sigma8 / const * pk(ik) 
      enddo
c
      write(6,*) 'const=', const
c
      end
c
c
c ******************************************************* c
c
      subroutine truncation_k(kmin, kmax)
c
c ******************************************************* c
c
      implicit none
c
      integer ik, ikk, ikmax, ik_max, ibox
      parameter(ikmax=3000)
      real*8 ak(ikmax), pk(ikmax)
      real*8 akk(ikmax), pkk(ikmax)
      real*8 kmin, kmax, Lbox, pi
      common /normalized_pk/ ak, pk, ik_max
      pi = 4.d0 * datan(1.d0)
c     -----------------------------------------------------
c
      do ik=1, ik_max
         akk(ik) = ak(ik)
         pkk(ik) = pk(ik)
      enddo
c
      ikk = 1
      do ik=1, ik_max
         if(akk(ik).ge.kmin .and. akk(ik).le.kmax) then
            ak(ikk) = akk(ik)
            pk(ikk) = pkk(ik)
            ikk = ikk + 1
         endif
      enddo
c
      ik_max = ikk -1
c
      write(6,*) 'ak(1)=', ak(1)
      write(6,*) 'ak(ik_max)=', ak(ik_max)
      write(6,*) 'ik_max=', ik_max
c
      end
c
c ******************************************************* c
c
      subroutine find_pk(kk, pklin)
c
c ******************************************************* c
c
      implicit none
      integer ik_max, ikmax 
      integer j, jmin, jmax
      parameter(ikmax=3000)
      real*8 ak(ikmax), pk(ikmax), kk, s, ds, pklin
      common /normalized_pk/ ak, pk, ik_max
c     -------------------------------------------
c
      call hunt(ak, ik_max, kk, j)
c
      jmin = j - 2
      jmax = j + 2
      if(jmin.lt.1) jmin = 1
      if(jmax.ge.ik_max) jmax = ik_max
c
      call polint(ak(jmin),pk(jmin),jmax-jmin+1,kk,s,ds)
      pklin = s
c      
      end
c
c ******************************************************* c
c
      subroutine calc_one_loop_pk(pk13, pk22)
c
c ******************************************************* c
c
      implicit none
      integer ik, ik_max, ikmax
      integer ix, ix_max, imu, imu_max
      parameter(ikmax=3000)
      parameter(ix_max=1000, imu_max=50)
      real*8  r_th, pi
      real*8  ak(ikmax), pk(ikmax), pk13(ikmax), pk22(ikmax) 
      real*8  kmin, kmax, xmin, xmax, mumin, mumax
      real*8  k, mu, fp13, fp22
      real*8  x1, x2, pklin1, pklin2, pklink
      real*8  mu1, mu2, ar1, ar2, pklin_ar1, pklin_ar2
      real*8  integrand1, integrand2
      common /normalized_pk/ ak, pk, ik_max
      pi = 4.d0 * datan(1.d0)
c     ---------------------------------------------------
c
      kmin = ak(1)
      kmax = ak(ik_max) 
c
      do 10 ik=1, ik_max
c
         k = ak(ik)
         call find_pk(k, pklink) 
c
         xmin = kmin / k
         xmax = kmax / k
c
         pk13(ik) = 0.0
         pk22(ik) = 0.0
         integrand1 = 0.0 
         integrand2 = 0.0 
c
         x1 = xmin 
         call find_pk(k*x1, pklin1)
c
         do 20 ix=1, ix_max-1
c
            x2 = xmin * (xmax/xmin) ** (dble(ix)/dble(ix_max-1))
            call find_pk(k*x2, pklin2)
c
ccc            write(6,*) ix
c
c ////// integration of P_13 ////// c
c
            pk13(ik) = pk13(ik) + 
     &           (fp13(x1) * pklin1 +  fp13(x2) * pklin2)* (x2 - x1)/2. 
c
c ////// integration of P_22 ////// c
c
            if(ix.eq.1) then
               integrand1 = 0.0
               mumin = max( -1.0, (1.+x1**2-xmax**2)/2./x1 ) 
               if(x1.le.0.5) mumax = min(1.0, (1.+x1**2-xmin**2)/2./x1)
               if(x1.gt.0.5) mumax = 0.5/x1
c
               mu1 = mumin 
               ar1 = k * dsqrt(1.+x1**2-2.*mu1*x1)
               call find_pk(ar1, pklin_ar1)
               do 25 imu=1, imu_max-1
                  mu2 = mumin + (mumax-mumin)*dble(imu)/dble(imu_max-1)
                  ar2 = k * dsqrt(1.+x1**2-2.*mu2*x1)
                  call find_pk(ar2, pklin_ar2)
                  integrand1 = integrand1 + 
     &                 ( pklin_ar1*fp22(mu1,x1)+pklin_ar2*fp22(mu2,x1) )
     &                 * (mu2 - mu1)/2.
                  mu1 = mu2
                  ar1 = ar2
                  pklin_ar1 = pklin_ar2
 25            continue
            endif
c
            mumin = max( -1.0, (1.+x2**2-xmax**2)/2./x2 )
            if(x2.le.0.5) mumax = min(1.0, (1.+x2**2-xmin**2)/2./x2)
            if(x2.gt.0.5) mumax = 0.5/x2
c
            integrand2 = 0.0
            mu1 = mumin 
            ar1 = k * dsqrt(1.+x2**2-2.*mu1*x2)
            call find_pk(ar1, pklin_ar1)
cccc
cccc            write(6,*) ar1, kmin, kmax, mumin, mumax, 
cccc     &           pklin_ar1
cccc
c
            do 30 imu=1, imu_max-1
               mu2 = mumin + (mumax-mumin)*dble(imu)/dble(imu_max-1)
               ar2 = k * dsqrt(1.+x2**2-2.*mu2*x2)
               call find_pk(ar2, pklin_ar2)
               integrand2 = integrand2 + 
     &              ( pklin_ar1*fp22(mu1,x2) + pklin_ar2*fp22(mu2,x2) )
     &              * (mu2 - mu1)/2.
               mu1 = mu2
               ar1 = ar2
               pklin_ar1 = pklin_ar2
 30         continue
            pk22(ik) = pk22(ik) + 
     &           (integrand1*pklin1 + integrand2*pklin2) * (x2 - x1)/2. 
            x1 = x2
            pklin1 = pklin2
            integrand1 = integrand2
 20   continue
c
      pk13(ik) = pklink * pk13(ik) * k**3 / (2.*pi)**2 
      pk22(ik) = 2. * pk22(ik) * k**3 / (2.*pi)**2 
c
      write(6,'(1p3e18.10)') ak(ik), pk13(ik), pk22(ik)
c      
 10   continue
c
      end
c
c ******************************************************* c
c
      function  fp13(x)
c
c ******************************************************* c
c
      implicit none
      real*8 fp13, x, s
c
      if(x.lt.0.) then
         write(6,*) ' fp13 negative !! emergency stop ' 
         stop
      elseif(x .lt. 5.d-3) then
         fp13 = -2./3. + 232./315.*x*x - 376./735.*x*x*x*x 
      elseif(x .gt. 5.d2) then
         s = 1.d0 / x
         fp13 = -122./315. + 8./105.*s*s -40./1323.*s*s*s*s
      elseif(x.ge.0.995 .and. x.le.1.005) then
         fp13 = ( -22. + 2.*(x-1.) -29.*(x-1.)**2 )/63.
      else 
         fp13 = ( 12./x/x - 158. + 100.*x*x -42.*x*x*x*x + 
     &        3./(x*x*x) * (x*x-1.)**3 * (7.*x*x + 2.) * 
     &        dlog(abs((x+1.)/(x-1.))) )/252.
      endif
c
      end
c
c ******************************************************* c
c
      function  fp22(mu,x)
c
c ******************************************************* c
c
      implicit none
      real*8 fp22, mu, x
c
      if(x.eq.1. .and. mu.eq.1.0) then
         write(6,*) ' fp22 pole !! integration stopped '
         stop
      endif
c
      fp22 = (3.*x + 7.*mu - 10.*mu*mu*x)/(1. + x*x - 2.*mu*x)/7.
      fp22 = fp22 * fp22 /2.
c
      end
c
c ******************************************************* c
c
      subroutine calc_correction(pk_corr)
c
c ******************************************************* c
c
      implicit none
      integer ik, ik_max, ikmax, ix, ix_max
      parameter(ikmax=3000)
      parameter(ix_max=3000)
      real*8  pi, x1, x2, pklin1, pklin2
      real*8  ak(ikmax), pk(ikmax), pk_corr(ikmax)
      real*8  k, kmin, kmax
      common /normalized_pk/ ak, pk, ik_max
      pi = 4.d0 * datan(1.d0)
c     ---------------------------------------------------
c
      kmin = ak(1)
      kmax = ak(ik_max) 
c
      do 10 ik=1, ik_max
c
         k = ak(ik)
         pk_corr(ik) = 0.0
         x1 = kmin 
         call find_pk(x1, pklin1)
c
         do 20 ix=1, ix_max-1
c
            x2 = kmin * (kmax/kmin) ** (dble(ix)/dble(ix_max-1))
            call find_pk(x2, pklin2)
c
            pk_corr(ik) = pk_corr(ik) + 
     &           (pklin1 +  pklin2)* (x2 - x1)/2. 
c
            x1 = x2
            pklin1 = pklin2
 20   continue
c
      pk_corr(ik) = pk_corr(ik) * k**2 / (6.d0 * pi**2)
c
      write(6,'(1p2e18.10)') ak(ik), pk_corr(ik)
c      
 10   continue
c
      end
c
c ************************************************************
c
      SUBROUTINE polint(xa,ya,n,x,y,dy)
c
c ************************************************************
c
      implicit none
      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
C
c ************************************************ c
c
      SUBROUTINE hunt(xx,n,x,jlo)
c
c ************************************************ c
c
      implicit none
      INTEGER jlo,n
      REAL*8 x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd
      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END
c
