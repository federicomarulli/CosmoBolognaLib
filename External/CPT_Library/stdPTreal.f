c ******************************************************* c
c                                                         c
c      (Standard) Perturbation Theory for                 c
c           Non-linver Evolution of Power Spectrum        c
c             using Gaussian quadratures                  c
c                                                         c  
c                                                         c  
c   The code provides the PT correction data, pk13 and    c
c   pk22, which are used for calculating non-linear       c
c   1-loop P(k). To compute 1-loop matter P(k) at         c
c   arbitrary redshift, you should also run the common    c
c   read file, read_pk.f                                  c  
c                                                         c  
c            Time-stamp: <2010-06-28 14:24:36 ataruya>    c
c ******************************************************* c
c
c     Note--. 
c
c     All the parameter 'ikmax' in this program must be 
c     the same. 
c
      program stdPTreal
c
      implicit none
c
      integer  ik, ikmax, ik_max, itype
      parameter(ikmax=3000)
      real*8   ak(ikmax), pk(ikmax), pk22(ikmax), pk13(ikmax)
      real*8   pi, sigma8, kmin, kmax
      character infile*50
      common /normalized_pk/ ak, pk, ik_max
      pi = 4.d0 * atan(1.d0)
c     -------------------------------------------------
c
      write(6,*)
      write(6,*) '**** Perturbation calculation of power spectrum ****'
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
c     /////// Save output data ///////
c
      open(10,file='stdPT.dat',status='unknown')
c
      do ik =1, ik_max
         write(10,'(1p4e18.10)') ak(ik), pk(ik), pk13(ik), pk22(ik)
      enddo
c
      write(6,*)
      write(6,*) ' Save output data, k, pk_lin, pk^(13), pk^(22), '
      write(6,*) '                                   to stdPT.dat '
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
      end
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
      integer ik, ik_max, ikmax, isub
      integer ix, ixmax13, ixmax22
      parameter(ikmax=3000)
      parameter(ixmax13=100, ixmax22=600)
      real*8  pi
      real*8  ak(ikmax), pk(ikmax), pk13(ikmax), pk22(ikmax)
      real*8  kmin, kmax, xmin, xmax, mumin, mumax
      real*8  k, pklink
      real*8  w13(ixmax13), x13(ixmax13), w22(ixmax22), x22(ixmax22)
      real*8  fp13, integ_fp22
      common /normalized_pk/ ak, pk, ik_max
      common /wave_number/  k, xmin, xmax
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
         pk13(ik)= 0.0
         pk22(ik)= 0.0

         xmin = kmin / k
         xmax = kmax / k
c
c     ////// Gauss-Legendre integration of p13(k), p22(k) //////  c
c
         call gauleg(log(xmin),log(xmax),x13,w13,ixmax13)
c
         do ix=1, ixmax13
            x13(ix)= dexp(x13(ix))
            pk13(ik) = pk13(ik) + w13(ix) * fp13(x13(ix))
         enddo
c
         if(k.lt.0.2) isub =200 
         if(k.ge.0.2) isub =0 
c
            call gauleg(log(xmin),log(xmax),x22,w22,ixmax22-isub)
c
            do ix=1, ixmax22-isub
               x22(ix)= dexp(x22(ix))
               pk22(ik) = pk22(ik) + w22(ix) * integ_fp22(x22(ix))
            enddo
c
c
         pk13(ik) = pk13(ik) * pklink * k**3 / (2.*pi)**2 
         pk22(ik) = 2.d0 * pk22(ik) * k**3 / (2.*pi)**2 
c
c
      write(6,'(1p4e18.10)') k, pklink, pk13(ik), pk22(ik)
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
      real*8 fp13, x, k, s, pklin
      real*8 xmin, xmax
      common /wave_number/  k, xmin, xmax
c
      fp13 = 0.0
c
      if(x.lt.0.) then
         write(6,*) ' warning: fp13 negative !! : x=',x
         fp13 = 0.0
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
      call find_pk(k*x, pklin)
c
      fp13 = fp13 * pklin * x
c
      end
c
c ******************************************************* c
c
      function integ_fp22(x)
c
c ******************************************************* c
c
      implicit none
      integer imu, imu_max
      parameter(imu_max=10)
      real*8  integ_fp22, xmin, xmax, mumin, mumax, fp22
      real*8  k, x, mu, wmu(imu_max), mmu(imu_max) 
      common /wave_number/  k, xmin, xmax
c
      integ_fp22 = 0.d0
c
      mumin = max(-1.0, (1.+x**2-xmax**2)/2./x)
      mumax = min( 1.0, (1.+x**2-xmin**2)/2./x)
c
      if(x.ge.0.5d0) mumax= 0.5d0/x
c
      call gauleg(mumin, mumax, mmu, wmu, imu_max)
c
      do imu=1, imu_max
         integ_fp22 = integ_fp22 + wmu(imu) * fp22(x, mmu(imu))
      enddo
c
      end
c
c ******************************************************* c
c
      function  fp22(x, mu)
c
c ******************************************************* c
c
      implicit none
      real*8 fp22
      real*8 mu, k, x, pklin1, pklin2  
      real*8 mumin, mumax, xmax, xmin 
      common /wave_number/  k, xmin, xmax
c
      fp22 = (3.*x + 7.*mu - 10.*mu*mu*x)/(1. + x*x - 2.*mu*x)/7.
c
      call find_pk(k*x, pklin1)
      call find_pk(k*sqrt(1.+x*x-2.*mu*x), pklin2)
c
      fp22 = fp22 * fp22 /2. * x * pklin1 * pklin2
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
c ************************************************************
c
      SUBROUTINE gauleg(x1,x2,x,w,n)
c
c ************************************************************
c
      INTEGER n
      REAL*8 x1,x2,x(n),w(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.d-14)
      INTEGER i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software z!0(0.
