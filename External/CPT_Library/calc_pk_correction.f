c ******************************************************* c
c                                                         c
c      Corrections to Scoccimarro (2004) formula          c
c      in redshift-space power spectrum                   c
c                                                         c  
c            Time-stamp: <2016-02-26 18:51:23 ataruya>    c
c                                                         c
c  inconsistency in common block /no_wiggle_param/ fixed  c
c ******************************************************* c
c
c     Note--. 
c
c     All the parameter 'ikmax' in this program must be 
c     the same. 
c
      program calc_pk_correction
c
      implicit none
c
      integer  ik, ikmax, ik_max, ik_max1, ik_max2
      integer  itype, ibox, isigmav
      parameter(ikmax=3000)
      real*8  ak(ikmax), pk(ikmax)
      real*8  akcorr(ikmax)
      real*8  pk0corr(ikmax), pk2corr(ikmax), pk4corr(ikmax)
      real*8  pk_B1(ikmax), pk_B2(ikmax), pk_B3(ikmax), pk_B4(ikmax)
      real*8  akEH(ikmax), pk0EH(ikmax), pk2EH(ikmax), pk4EH(ikmax)
      real*8   pi, Lbox, zred, ss, sigmav, growth, f
      character infile*50
      common /pk_data/ ak, pk, ik_max
      common /pkcorr/ akcorr, pk0corr, pk2corr, pk4corr, 
     &     pk_B1, pk_B2, pk_B3, pk_B4, ik_max1
      common /pkEH/ akEH, pk0EH, pk2EH, pk4EH, ik_max2
      pi = 4.d0 * atan(1.d0)
c     -------------------------------------------------
c
      write(6,*)
      write(6,*) '*** Corrections to Scoccimarro (2004) formula ***'
      write(6,*) '***     in redshift-space power spectrum      ***'
      write(6,*)
      write(6,*)
      write(6,*) 'Output redshift ?'
      read(5,*) zred
      write(6,*) 'including finite-volume effect? :type (0 or 1, L_box)'
      write(6,*) ' (e.g., (0, 300) for boxsize L_box=300Mpc/h ) '
      write(6,*) ' (      (1, ???) for ignoring finite-volume effect) '
      read(5,*) ibox, Lbox
c
c     /////// Load (linear) matter power spectrum ///////
c
      call load_matterpower_data
c
      write(6,*) ' loading (linear) matter power spectrum, done '
c
c     //////// Set cosmological parameters ////////
c
      call set_cosmological_params
c
c     /////// Sigma8 normalization ///////
c
      call normalization_trapez
c
      write(6,*) ' normalization done '
c
c     //////// Truncation of low-k, high-k  in linear P(k) ////////
c
      write(6,*)
      write(6,'(A,1p1e12.6)') 'f(z)=',f(zred)
      write(6,*)
c
      call calc_sigmav(ss)
      sigmav = ss
      write(6,'(A,1p1e12.6)') 'sigma_v=',
     &     ss * growth(zred)
      write(6,*)
      write(6,*) 'Use this value ? y[0], n[1]'
      read(5,*) isigmav
      if(isigmav.eq.1) then
         write(6,*) 'Input sigma_v'
         read(5,*) sigmav
         sigmav = sigmav / growth(zred)
      endif
c
      call truncation_k(ibox, Lbox)
c
c     /////// Correction to Scoccimarro (2004) formula ///////
c
      call calc_correction
c
      write(6,*) ' 1-loop P(k) done '
c
c     /////// Summing up all contributions ///////
c
      call calc_pkred(zred,sigmav)
c
      write(6,*) ' summing up all contributions done '
c
c     /////// Save output data ///////
c
      open(10,file='corr_pkred.dat',status='unknown')
      open(11,file='pkcorr_red.dat',status='unknown')
c
      do ik =1, ik_max
         write(10,'(1p8e18.10)') ak(ik), pk(ik), 
     &        pk0EH(ik), pk0corr(ik), pk2EH(ik), 
     &        pk2corr(ik), pk4EH(ik), pk4corr(ik) 
         write(11,'(1p5e18.10)') ak(ik), pk_B1(ik), pk_B2(ik), 
     &        pk_B3(ik), pk_B4(ik)
      enddo
c
      close(10)
      close(11)
c
      end
c
c ******************************************************* c
c
      subroutine load_matterpower_data
c
c ******************************************************* c
c
c     input file is assumed to be the matter power spectrum data 
c     created by CAMB code. 
c
      implicit none
c
      integer ik, ikk, ikmax, ik_max
      parameter(ikmax=3000)
      real*8  ak(ikmax), pk(ikmax), pi
      common /pk_data/ ak, pk, ik_max
      character infile*50
      pi = 4.d0 * datan(1.d0)
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
      write(6,*)
c
      close(9)
c
      end
c
c ******************************************************* c
c
      subroutine set_cosmological_params
c
c ******************************************************* c
c
c     Set cosmological parameters for linear growth factor,
c     growth-rate parameter, and sigma8 normalization 
c
c     Note.--
c     assume a flat cosmological model (Omega_v = 1-Omega_m)
c
      implicit none
      integer iparams
      real*8 h, Tcmb, n_s, sigma8, Omega_m, Omega_b, Omega_v, w_de
      common /cosmological_param/ Omega_m, Omega_v, w_de, sigma8
      common /no_wiggle_param/ Omega_b, h, Tcmb, n_s
c     ---------------------------------------------------
c
      write(6,*)
      write(6,*) '*** Set cosmological parameters ***' 
      write(6,*)
c
c     Default cosmological parameters (you can change them appropriately)
c
      h = 0.701d0
      Tcmb = 2.726d0
      n_s = 0.96d0
      sigma8 = 0.817d0
      Omega_m = 0.279d0
      Omega_b = 0.165 * Omega_m
      Omega_v = 1.d0 - Omega_m
      w_de = -1.d0
c
 5    write(6,*) '[1] h        =',h
      write(6,*) '[2] T_cmb[K] =',Tcmb
      write(6,*) '[3] n_s      =',n_s
      write(6,*) '[4] sigma_8  =',sigma8
      write(6,*) '[5] Omega_m  =',Omega_m
      write(6,*) '[6] Omega_b  =',Omega_b
      write(6,*) '[7] w_de     =',w_de
      write(6,*) 
      write(6,*) 'Note--. spatial curvature is forced to be flat, '
      write(6,*) '                   i.e., Omega_v = 1 - Omega_m  '
      write(6,*)
      write(6,*) 'change cosmological parameter? [1-7] or n[0]'
      read(5,*)  iparams
      if(iparams.eq.0) goto 8
      if(iparams.eq.1) then
         write(6,*) 'type h'
         read(5,*) h
      elseif(iparams.eq.2) then
         write(6,*) 'type Tcmb'
         read(5,*) Tcmb
      elseif(iparams.eq.3) then
         write(6,*) 'type n_s'
         read(5,*) n_s
      elseif(iparams.eq.4) then
         write(6,*) 'type sigma8'
         read(5,*) sigma8
      elseif(iparams.eq.5) then
         write(6,*) 'type Omega_m'
         read(5,*) Omega_m
      elseif(iparams.eq.6) then
         write(6,*) 'type Omega_b'
         read(5,*) Omega_b
      elseif(iparams.eq.7) then
         write(6,*) 'type w'
         read(5,*) w_de
      else 
         stop
      endif
      goto 5
c
 8    continue
c
      end
c
c ******************************************************* c
c
      subroutine normalization_trapez
c
c ******************************************************* c
c
      implicit none
      integer ik, ik_max, ikmax
      parameter(ikmax=3000)
      real*8  Omega_m, Omega_v, w_de, sigma8
      real*8  ak(ikmax), pk(ikmax)
      real*8  r_th, const, x
      real*8  W_TH, sigma_a, sigma_b, pi
      common /pk_data/ ak, pk, ik_max
      common /cosmological_param/ Omega_m, Omega_v, w_de, sigma8
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
      subroutine truncation_k(ibox, Lbox)
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
      common /pk_data/ ak, pk, ik_max
      pi = 4.d0 * datan(1.d0)
c     -----------------------------------------------------
c
      kmin = 1.d-5   ! default value
c
      if(ibox.eq.0) kmin = 2.d0 * pi / Lbox
c
      kmax = 600.d0    ! default value
c
      do ik=1, ik_max
         akk(ik) = ak(ik)
         pkk(ik) = pk(ik)
      enddo
c
      ikk = 1
      do ik=1, ik_max
         if(akk(ik).ge.kmin .and. akk(ik).le.kmax
cc     &        .and. mod(ik,2).eq.0 ) then
     &        ) then
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
      write(6,*)
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
      common /pk_data/ ak, pk, ik_max
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
      subroutine calc_correction
c
c ******************************************************* c
c
      implicit none
      integer ik, ik_max, ikmax, isub
      integer ix, ixmax
      parameter(ikmax=3000)
      parameter(ixmax=600)
      real*8  pi
      real*8  ak(ikmax), pk(ikmax)
      real*8  pk_B111(ikmax), pk_B112(ikmax), pk_B121(ikmax)
      real*8  pk_B122(ikmax), pk_B211(ikmax), pk_B212(ikmax)
      real*8  pk_B221(ikmax), pk_B222(ikmax), pk_B312(ikmax)
      real*8  pk_B321(ikmax), pk_B322(ikmax), pk_B422(ikmax)
      real*8  kmin, kmax, xmin, xmax, mumin, mumax
      real*8  k, ww(ixmax), xx(ixmax), integ_fp
      common /pk_data/ ak, pk, ik_max
      common /wave_number/  k, xmin, xmax
      common /corr_pk/ pk_B111, pk_B112, pk_B121, pk_B122, pk_B211, 
     &     pk_B212, pk_B221, pk_B222, pk_B312, pk_B321, pk_B322, pk_B422
      pi = 4.d0 * datan(1.d0)
c     ---------------------------------------------------
c
      kmin = ak(1)
      kmax = ak(ik_max) 
c
      do 10 ik=1, ik_max
c
         k = ak(ik)
         pk_B111(ik)= 0.d0 
         pk_B112(ik)=0.d0
         pk_B121(ik)=0.d0
         pk_B122(ik)=0.d0
         pk_B211(ik)=0.d0
         pk_B212(ik)=0.d0
         pk_B221(ik)=0.d0
         pk_B222(ik)=0.d0
         pk_B312(ik)=0.d0
         pk_B321(ik)=0.d0
         pk_B322(ik)=0.d0
         pk_B422(ik)=0.d0
c
         xmin = kmin / k
         xmax = kmax / k
c
c     ////// Gauss-Legendre integration //////  c
c
         if(k.lt.0.2) isub =200 
         if(k.ge.0.2) isub =0 
c
            call gauleg(log(xmin),log(xmax),xx,ww,ixmax-isub)
c
            do ix=1, ixmax-isub
               xx(ix)= dexp(xx(ix))
               pk_B111(ik) = pk_B111(ik)+ww(ix)*integ_fp(1,xx(ix))
               pk_B112(ik) = pk_B112(ik)+ww(ix)*integ_fp(2,xx(ix))
               pk_B121(ik) = pk_B121(ik)+ww(ix)*integ_fp(3,xx(ix))
               pk_B122(ik) = pk_B122(ik)+ww(ix)*integ_fp(4,xx(ix))
               pk_B211(ik) = pk_B211(ik)+ww(ix)*integ_fp(5,xx(ix))
               pk_B212(ik) = pk_B212(ik)+ww(ix)*integ_fp(6,xx(ix))
               pk_B221(ik) = pk_B221(ik)+ww(ix)*integ_fp(7,xx(ix))
               pk_B222(ik) = pk_B222(ik)+ww(ix)*integ_fp(8,xx(ix))
               pk_B312(ik) = pk_B312(ik)+ww(ix)*integ_fp(9,xx(ix))
               pk_B321(ik) = pk_B321(ik)+ww(ix)*integ_fp(10,xx(ix))
               pk_B322(ik) = pk_B322(ik)+ww(ix)*integ_fp(11,xx(ix))
               pk_B422(ik) = pk_B422(ik)+ww(ix)*integ_fp(12,xx(ix))
            enddo
c
         pk_B111(ik) = 2.d0 * pk_B111(ik) * k**3 / (2.*pi)**2
         pk_B112(ik) = - 2.d0 * pk_B112(ik) * k**3 / (2.*pi)**2
         pk_B121(ik) = - 2.d0 * pk_B121(ik) * k**3 / (2.*pi)**2
         pk_B122(ik) = 2.d0 * pk_B122(ik) * k**3 / (2.*pi)**2
         pk_B211(ik) = 2.d0 * pk_B211(ik) * k**3 / (2.*pi)**2
         pk_B212(ik) = - 2.d0 * pk_B212(ik) * k**3 / (2.*pi)**2
         pk_B221(ik) = - 2.d0 * pk_B221(ik) * k**3 / (2.*pi)**2
         pk_B222(ik) = 2.d0 * pk_B222(ik) * k**3 / (2.*pi)**2
         pk_B312(ik) = - 2.d0 * pk_B312(ik) * k**3 / (2.*pi)**2
         pk_B321(ik) = - 2.d0 * pk_B321(ik) * k**3 / (2.*pi)**2
         pk_B322(ik) = 2.d0 * pk_B322(ik) * k**3 / (2.*pi)**2
         pk_B422(ik) = 2.d0 * pk_B422(ik) * k**3 / (2.*pi)**2
c
      write(6,'(i3,1p3e18.10)') ik,k,pk_B111(ik),pk_B112(ik)
c
 10   continue
c
      end
c
c ******************************************************* c
c
      function integ_fp(ip, x)
c
c ******************************************************* c
c
      implicit none
      integer ip, imu, imu_max
      parameter(imu_max=10)
      real*8  integ_fp, xmin, xmax, mumin, mumax, fp
      real*8  k, x, mu, wmu(imu_max), mmu(imu_max)
      common /wave_number/  k, xmin, xmax
c     ------------------------------------------  c
c
      integ_fp = 0.d0
c
      mumin = max(-1.0, (1.+x**2-xmax**2)/2./x)
      mumax = min( 1.0, (1.+x**2-xmin**2)/2./x)
c
      if(x.ge.0.5d0) mumax= 0.5d0/x
c
      call gauleg(mumin, mumax, mmu, wmu, imu_max)
c
      do imu=1, imu_max
         integ_fp = integ_fp + wmu(imu) * fp(ip, x, mmu(imu))
      enddo
c
      end
c
c ******************************************************* c
c
      function  fp(ip, x, mu)
c
c ******************************************************* c
c
c     ip=1 for kernel of pk_B111
c     ip=2 for kernel of pk_B112
c     ip=3 for kernel of pk_B121
c     ip=4 for kernel of pk_B122
c     ip=5 for kernel of pk_B211
c     ip=6 for kernel of pk_B212
c     ip=7 for kernel of pk_B221
c     ip=8 for kernel of pk_B222
c     ip=9 for kernel of pk_B312
c     ip=10 for kernel of pk_B321
c     ip=11 for kernel of pk_B322
c     ip=12 for kernel of pk_B422
c
      implicit none
      integer ip
      real*8 fp
      real*8 mu, k, x, pklin1, pklin2
      real*8 mumin, mumax, xmax, xmin
      common /wave_number/  k, xmin, xmax
c     --------------------------------------- c
c
      if(ip.eq.1) then
         fp = x**2 * (mu*mu-1.) / 2.
      elseif(ip.eq.2) then
         fp = 3.*x**2 * (mu*mu-1.)**2 / 8.
      elseif(ip.eq.3) then
         fp = 3.*x**4 * (mu*mu-1.)**2 / (1.+x*x-2.*mu*x) / 8.
      elseif(ip.eq.4) then
         fp = 5.*x**4 * (mu*mu-1.)**3 / (1.+x*x-2.*mu*x) / 16.
      elseif(ip.eq.5) then
         fp = x * (x+2.*mu-3.*x*mu*mu) / 2.
      elseif(ip.eq.6) then
         fp = - 3.*x * (mu*mu-1.) * (-x-2.*mu+5.*x*mu*mu) / 4.
      elseif(ip.eq.7) then
         fp = 3.*x**2 * (mu*mu-1.) * (-2.+x*x+6.*x*mu-5.*x*x*mu*mu) 
     &        / (1.+x*x-2.*mu*x) / 4.
      elseif(ip.eq.8) then
         fp = - 3.*x**2 * (mu*mu-1.)**2 
     &        * (6.-5.*x*x-30.*x*mu+35.*x*x*mu*mu) 
     &        / (1.+x*x-2.*mu*x) / 16.
      elseif(ip.eq.9) then
         fp = x * (4.*mu*(3.-5.*mu*mu) + x*(3.-30.*mu*mu+35.*mu**4) )
     &        / 8.
      elseif(ip.eq.10) then
         fp = x * (-8.*mu + x*(-12.+36.*mu*mu+12.*x*mu*(3.-5.*mu*mu)+
     &        x**2*(3.-30.*mu*mu+35.*mu**4) ) ) / (1.+x*x-2.*mu*x) / 8.
      elseif(ip.eq.11) then
         fp = 3.*x * (mu*mu-1.) * (-8.*mu + x*(-12.+60.*mu*mu+
     &        20.*x*mu*(3.-7.*mu*mu)+5.*x*x*(1.-14.*mu*mu+21.*mu**4)) )
     &        / (1.+x*x-2.*mu*x) / 16.
      elseif(ip.eq.12) then
         fp = x * (8.*mu*(-3.+5.*mu*mu) - 6.*x*(3.-30.*mu*mu+35.*mu**4) 
     &        + 6.*x*x*mu*(15.-70.*mu*mu+63*mu**4) + x**3*(5.-21.*mu*mu*
     &        (5.-15.*mu*mu+11.*mu**4)) ) / (1.+x*x-2.*mu*x) / 16.
      endif
c
      call find_pk(k*x, pklin1)
      call find_pk(k*sqrt(1.+x*x-2.*mu*x), pklin2)
c
      fp = fp * x * pklin1 * pklin2 / (1.+x*x-2.*mu*x)
c
      end
c
c ******************************************************* c
c
      subroutine calc_pkred(zred, sigmav)
c
c ******************************************************* c
c
c     Summing up all contributions to redshift P(k) in PT
c     and calculating monopole, quadrupole and hexadecapole
c     spectra
c
      implicit none
      integer ik, ikmax, ik_max, ik_max1, ik_max2
      parameter(ikmax=3000)
      real*8  zred, growth, f, dd2, ff
      real*8  Omega_m, Omega_v, w_de, sigma8
      real*8  Omega_b, h, Tcmb, n_s
      real*8  ak(ikmax), pk(ikmax), pk_EH(ikmax)
      real*8  pk_B111(ikmax), pk_B112(ikmax), pk_B121(ikmax)
      real*8  pk_B122(ikmax), pk_B211(ikmax), pk_B212(ikmax)
      real*8  pk_B221(ikmax), pk_B222(ikmax), pk_B312(ikmax)
      real*8  pk_B321(ikmax), pk_B322(ikmax), pk_B422(ikmax)
      real*8  akcorr(ikmax)
      real*8  pk0corr(ikmax), pk2corr(ikmax), pk4corr(ikmax)
      real*8  pk_B1(ikmax), pk_B2(ikmax), pk_B3(ikmax), pk_B4(ikmax)
      real*8  akEH(ikmax), pk0EH(ikmax), pk2EH(ikmax), pk4EH(ikmax)
      real*8  fact, alpha, sigmav
      common /alpha_param/ alpha
      common /cosmological_param/ Omega_m, Omega_v, w_de, sigma8
      common /no_wiggle_param/ Omega_b, h, Tcmb, n_s
      common /pk_data/ ak, pk, ik_max
      common /corr_pk/ pk_B111, pk_B112, pk_B121, pk_B122, pk_B211, 
     &     pk_B212, pk_B221, pk_B222, pk_B312, pk_B321, pk_B322, pk_B422
      common /pkcorr/ akcorr, pk0corr, pk2corr, pk4corr, 
     &     pk_B1, pk_B2, pk_B3, pk_B4, ik_max1
      common /pkEH/ akEH, pk0EH, pk2EH, pk4EH, ik_max2
c     -------------------------------------------------------- c
      ik_max1 = ik_max
      ik_max2 = ik_max
c
c     growth factor and its logarithmic derivative
c     (assuming flat universe)
      dd2 = growth(zred)**2
      ff = f(zred)
c
      write(6,*) 'd, f, sigmav=', dd2, ff, sigmav*dsqrt(dd2)
c
c     ////// corrections to power spectrum  ////// c
c
ccc      open(11,file='fact.dat',status='unknown')
      open(11,file='pkstd_corr_tree.dat',status='unknown')
c
      do ik=1, ik_max1
c
         alpha = (ak(ik)*ff*sigmav)**2 * dd2
c
         akcorr(ik) = ak(ik)
c
         pk_B1(ik) = ff**2*pk_B111(ik) + ff**3*pk_B112(ik) + 
     &        ff**3*pk_B121(ik) + ff**4*pk_B122(ik)
         pk_B2(ik) = ff**2*pk_B211(ik) + ff**3*pk_B212(ik) + 
     &        ff**3*pk_B221(ik) + ff**4*pk_B222(ik)
         pk_B3(ik) = ff**3*pk_B312(ik) + ff**3*pk_B321(ik) +
     &        ff**4*pk_B322(ik)
         pk_B4(ik) = ff**4*pk_B422(ik)
c
         pk_B1(ik) = pk_B1(ik) * dd2*dd2
         pk_B2(ik) = pk_B2(ik) * dd2*dd2
         pk_B3(ik) = pk_B3(ik) * dd2*dd2
         pk_B4(ik) = pk_B4(ik) * dd2*dd2
c
         pk0corr(ik) = fact(1,0) * pk_B1(ik) 
     &        + fact(2,0) * pk_B2(ik) + fact(3,0) * pk_B3(ik)
     &        + fact(4,0) * pk_B4(ik)
c
         pk2corr(ik) = fact(1,2) * pk_B1(ik) 
     &        + fact(2,2) * pk_B2(ik) + fact(3,2) * pk_B3(ik)
     &        + fact(4,2) * pk_B4(ik)
c
         pk4corr(ik) = fact(1,4) * pk_B1(ik)
     &        + fact(2,4) * pk_B2(ik) + fact(3,4) * pk_B3(ik)
     &        + fact(4,4) * pk_B4(ik)
c
ccc      write(11,'(1p13e18.10)') alpha, fact(1,0), fact(2,0), fact(3,0), 
ccc     &        fact(4,0), fact(1,2), fact(2,2), fact(3,2), fact(4,2),
cc     &        fact(1,4), fact(2,4), fact(3,4), fact(4,4)
         write(11,'(1p10e18.10)') akcorr(ik), 
     &        pk_B111(ik) * dd2*dd2,
     &        ( pk_B112(ik) + pk_B121(ik) ) * dd2*dd2,
     &        pk_B122(ik) * dd2*dd2,
     &        pk_B211(ik) * dd2*dd2,
     &        ( pk_B212(ik) + pk_B221(ik) ) * dd2*dd2,
     &        pk_B222(ik) * dd2*dd2,
     &        ( pk_B312(ik) + pk_B321(ik) ) * dd2*dd2,
     &        pk_B322(ik) * dd2*dd2,
     &        pk_B422(ik) * dd2*dd2
c
      enddo
c
      close(11)
      write(6,*) ' Output file: pkstd_corr_tree.dat '
c
c     ////// no-wiggle linear power spectrum  ////// c
c
      call no_wiggle_Pk(sigma8, ik_max, ak, pk_EH)
c
      ff = dabs(ff)
c
      do ik=1, ik_max2
         akEH(ik) = ak(ik)
         pk(ik) = pk(ik) * dd2 
         pk0EH(ik) = (1.d0+2.d0/3.d0*ff+1.d0/5.d0*ff*ff)* pk_EH(ik)*dd2
         pk2EH(ik) = ( 4.d0/3.d0*ff + 4.d0/7.d0*ff*ff )* pk_EH(ik)*dd2
         pk4EH(ik) = 8.d0/35.d0*ff*ff * pk_EH(ik)*dd2
      enddo
c
      end
c
c ******************************************************* c
c
      subroutine calc_sigmav(sigmav)
c
c ******************************************************* c
c
      implicit none
      integer ikmax, ik, ik_max
      parameter(ikmax=3000)
      real*8  ak(ikmax), pk(ikmax), sigmav, pi
      common /pk_data/ ak, pk, ik_max
c     -------------------------------------------------
      pi = 4.d0 * datan(1.d0)
c
      sigmav = 0.d0
c
      do ik=1, ik_max-1
         sigmav = sigmav + (pk(ik+1) + pk(ik)) * 
     &        (ak(ik+1)-ak(ik)) /2.d0
      enddo
c
      sigmav = dsqrt(sigmav /(2*pi*pi) /3.d0)
c
      end
c
c ************************************************ c
c
      function fact(n, l)
c
c ************************************************ c
c
c     (2l+1)/2 * integ dmu  mu^(2n) * exp(-alpha*mu^2) * P_l(mu)
c
      implicit none
      integer n, l
      real*8 fact, nn, alpha, gamhalf, gammp
      common /alpha_param/ alpha
c     ---------------------------- c
      nn = dble(n)
c
      if(alpha.gt.0.05) then
c
         if(l.eq.0) then
            fact = gamhalf(n) * gammp(0.5+nn,alpha)
            fact = fact / alpha**(nn+0.5) / 4.d0
         elseif(l.eq.2) then
            fact = alpha * gamhalf(n) * gammp(0.5+nn,alpha) 
     &           - 3.d0 * gamhalf(n+1) * gammp(1.5+nn,alpha) 
            fact = fact / alpha**(nn+1.5) * (-5.d0/8.d0)
         elseif(l.eq.4) then
            fact = 12.*gamhalf(n)*gammp(0.5+nn,alpha)/alpha**(n+0.5)  
     &           -120.*gamhalf(n+1)*gammp(1.5+nn,alpha)/alpha**(n+1.5)
     &           +140.*gamhalf(n+2)*gammp(2.5+nn,alpha)/alpha**(n+2.5)
            fact = fact * 9./128.
         endif
c
      else
c
         if(l.eq.0) then
            fact = 1./(2.+4.*nn) - alpha/(6.+4.*nn) 
     &           + alpha**2/(20.+8.*nn) 
         elseif(l.eq.2) then
            fact = nn/(3.+8.*nn+4.*nn**2) 
     &           - (nn+1.)*alpha/(15.+16.*nn+4.*nn**2) 
     &           + (nn+2.)*alpha**2/(70.+48.*nn+8.*nn**2)
            fact = fact * 5.d0
         elseif(l.eq.4) then
            fact = dble(n*(n-1))/dble(15+46*n+36*n**2+8*n**3) 
     &           - dble(n*(n+1))/dble(105+142*n+60*n**2+8*n**3)*alpha
     &           + dble((n+1)*(n+2))/dble(315+286*n+84*n**2+8*n**3)
     &           *alpha**2/2.d0
            fact = fact * 18.d0
         endif
c         
      endif
c
      fact = fact * (1.d0 + (-1.d0)**(2.*nn)) 
c
      end
c
c ************************************************ c
c
      function gamhalf(n)
c
c ************************************************ c
c
c     Gamma(n+1/2) up to n=6
c
      integer n
      real*8 gamhalf, pi
      pi = 4.d0 *datan(1.d0)
c
      if(n.eq.0) gamhalf = 1.d0
      if(n.eq.1) gamhalf = 0.5d0
      if(n.eq.2) gamhalf = 0.75d0
      if(n.eq.3) gamhalf = 1.875d0
      if(n.eq.4) gamhalf = 6.5625d0
      if(n.eq.5) gamhalf = 29.53125d0
      if(n.eq.6) gamhalf = 162.421875d0
      if(n.eq.7) gamhalf = 1055.7421875d0
      if(n.eq.8) gamhalf = 7918.06640625d0
c
      gamhalf = gamhalf * dsqrt(pi)
c
      end
c
c ************************************************ c
c
      subroutine no_wiggle_Pk(sigma8, ik_max, ak, pk_EH)
c
c ************************************************ c
c
      implicit none
      integer ik, ikmax, ik_max
      parameter(ikmax=3000)
      real*8  ak(ikmax), pk_EH(ikmax)
      real*8  sigma8
      real*8  const, r_th, ks, kmin, kmax, s1, s2
      real*8  Pk_lin_EH
      external var, midpnt, midinf
      common /R_tophat/ r_th
c     ----------------------------------
c
c ///// normalization by sigma8 ///// c
c
      r_th = 8.d0
      ks = 1.d0 / r_th
      kmin = 0.00001 * ks
      kmax = 1.d3
c
      call qromo(var, kmin, ks, s1, midpnt)
      call qromo(var, ks, kmax, s2, midinf)
c
      const = sigma8**2 / (s1 + s2)
c
      do ik = 1, ik_max
         pk_EH(ik) = const * pk_lin_EH(ak(ik))
      enddo
c
      end
c ************************************************ c
c
      function var(k) 
c
c ************************************************ c
c
      implicit none
      real*8  var, k, x, w_th, Pk_lin_EH
      real*8  pi, r_th
      common /R_tophat/ r_th
c     -----------------------------------
c
      pi = 4.d0 * datan(1.d0)
c
      x = k * r_th
      w_th = 3.d0 * (sin(x)-x*cos(x)) / x**3 
      var = k**2 * w_th**2 * Pk_lin_EH(k) / (2.d0*pi**2)
c
      end
c
c ************************************************ c
c
      function Pk_lin_EH(k)
c
c ************************************************ c
c
c     compute un-normalized linear P(k) 
c     based on eq.(29) of Eisenstein & Hu (1998)
c     (no-wiggle approximation)  
c
      implicit none
      real*8 k, ss, alpha_gam, theta_cmb
      real*8 gamma_eff, q, L0, C0, T_EH, Pk_lin_EH
      real*8 omegab, omega0, h, Tcmb, n_s
!! AT 2014/5/16
      real*8 omega_v, w_de, sigma8
      common /cosmological_param/ omega0, omega_v, w_de, sigma8
      common /no_wiggle_param/ omegab, h, Tcmb, n_s
!! AT 2014/5/16
c     -----------------------------------------
c
c ///// fitting formula for no-wiggle P(k) (Eq.[29] of EH98)
c
      ss = 44.5 * h * dlog( 9.83 / (omega0*h*h) ) / 
     &     dsqrt( 1.d0 + 10.d0 * (omegab*h*h)**0.75 )
      alpha_gam = 1.d0 
     &     - 0.328 * dlog( 431. * omega0*h*h ) * omegab/omega0
     &     + 0.38 * dlog( 22.3 * omega0*h*h ) * (omegab/omega0)**2
      theta_cmb = Tcmb / 2.70 
      gamma_eff = omega0 * h * 
     &     ( alpha_gam + (1.d0 - alpha_gam) / (1.d0 + (0.43*k*ss)**4) )
c
      q = k * theta_cmb**2 / gamma_eff
      L0 = dlog( 2.d0 * dexp(1.d0) + 1.8 * q ) 
      C0 = 14.2 + 731.d0 / ( 1.d0 + 62.5 * q )
c
      T_EH = L0 / (L0 + C0*q*q )
c
      Pk_lin_EH = k ** n_s * T_EH**2 
c
      end
c
c ******************************************************* c
c
      function growth(zred)
c
c ******************************************************* c
c
c     Linear growth factor for flat cosmology
c
      implicit none
      real*8 growth, zred, Omega_m, Omega_v, w_de, sigma8, a, b, c
      real*8 zred1, zi, zz, gz, g0
      common /cosmological_param/ Omega_m, Omega_v, w_de, sigma8
c     --------------------------------------------
c
      a = -1.d0 / (3.d0 * w_de)
      b = (w_de - 1.d0)/ (2.d0 * w_de)
      c = 1.d0 - 5.d0 / (6.d0 * w_de)
c
      zred1 = 1.d0 + zred
      zi = - Omega_v / Omega_m
      zz = zi * zred1**(3.d0*w_de) 
c
      call HYGFX(a,b,c,zi, g0)  
      call HYGFX(a,b,c,zz, gz)  
c
      growth = (gz/g0) / zred1
c
      end
c
c ******************************************************* c
c
      function f(zred)
c
c ******************************************************* c
c
c     d lnD_+(z) / d ln a,   as function of redshift 
c
      implicit none
c
      real*8  f, zred, zred1, Omega_m, Omega_v, w_de, sigma8
      real*8  a, b, c, zi, zz, g1, g2
      common /cosmological_param/ Omega_m, Omega_v, w_de, sigma8
c     ---------------------------------------------------
c
      zred1 = 1.d0 + zred
      zi = - Omega_v / Omega_m
      zz = zi * zred1**(3.d0*w_de) 
c
      a = 1.d0 - 1.d0 / (3.d0 * w_de)
      b = 1.5d0 - 1.d0 / (2.d0 * w_de)
      c = 2.d0 - 5.d0 / (6.d0 * w_de)
      call HYGFX(a,b,c,zz, g1)  
c
      a = - 1.d0 / (3.d0 * w_de)
      b = (-1.d0 + w_de) / (2.d0 * w_de)
      c = 1.d0 - 5.d0 / (6.d0 * w_de)
      call HYGFX(a,b,c,zz, g2)  
c
      f = 1.d0 + 3.d0 * (-1.d0 + w_de) / (6.d0*w_de -5.d0) 
     &     * zz * g1 / g2
c
      end
c
c ******************************************************* c
c
        SUBROUTINE HYGFX(A,B,C,X,HF)
c
c ******************************************************* c
C
C       ====================================================
C       Purpose: Compute hypergeometric function F(a,b,c,x)
C       Input :  a --- Parameter
C                b --- Parameter
C                c --- Parameter, c <> 0,-1,-2,...
C                x --- Argument   ( x < 1 )
C       Output:  HF --- F(a,b,c,x)
C       Routines called:
C            (1) GAMMAX for computing gamma function
C            (2) PSI for computing psi function
C       ====================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        LOGICAL L0,L1,L2,L3,L4,L5
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        L0=C.EQ.INT(C).AND.C.LT.0.0
        L1=1.0D0-X.LT.1.0D-15.AND.C-A-B.LE.0.0
        L2=A.EQ.INT(A).AND.A.LT.0.0
        L3=B.EQ.INT(B).AND.B.LT.0.0
        L4=C-A.EQ.INT(C-A).AND.C-A.LE.0.0
        L5=C-B.EQ.INT(C-B).AND.C-B.LE.0.0
        IF (L0.OR.L1) THEN
           WRITE(*,*)'The hypergeometric series is divergent'
           RETURN
        ENDIF
        EPS=1.0D-15
        IF (X.GT.0.95) EPS=1.0D-8
        IF (X.EQ.0.0.OR.A.EQ.0.0.OR.B.EQ.0.0) THEN
           HF=1.0D0
           RETURN
        ELSE IF (1.0D0-X.EQ.EPS.AND.C-A-B.GT.0.0) THEN
           CALL GAMMAX(C,GC)
           CALL GAMMAX(C-A-B,GCAB)
           CALL GAMMAX(C-A,GCA)
           CALL GAMMAX(C-B,GCB)
           HF=GC*GCAB/(GCA*GCB)
           RETURN
        ELSE IF (1.0D0+X.LE.EPS.AND.DABS(C-A+B-1.0).LE.EPS) THEN
           G0=DSQRT(PI)*2.0D0**(-A)
           CALL GAMMAX(C,G1)
           CALL GAMMAX(1.0D0+A/2.0-B,G2)
           CALL GAMMAX(0.5D0+0.5*A,G3)
           HF=G0*G1/(G2*G3)
           RETURN
        ELSE IF (L2.OR.L3) THEN
           IF (L2) NM=INT(ABS(A))
           IF (L3) NM=INT(ABS(B))
           HF=1.0D0
           R=1.0D0
           DO 10 K=1,NM
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X
10            HF=HF+R
           RETURN
        ELSE IF (L4.OR.L5) THEN
           IF (L4) NM=INT(ABS(C-A))
           IF (L5) NM=INT(ABS(C-B))
           HF=1.0D0
           R=1.0D0
           DO 15 K=1,NM
              R=R*(C-A+K-1.0D0)*(C-B+K-1.0D0)/(K*(C+K-1.0D0))*X
15            HF=HF+R
           HF=(1.0D0-X)**(C-A-B)*HF
           RETURN
        ENDIF
        AA=A
        BB=B
        X1=X
        IF (X.LT.0.0D0) THEN
           X=X/(X-1.0D0)
           IF (C.GT.A.AND.B.LT.A.AND.B.GT.0.0) THEN
              A=BB
              B=AA
           ENDIF
           B=C-B
        ENDIF
        IF (X.GE.0.75D0) THEN
           GM=0.0D0
           IF (DABS(C-A-B-INT(C-A-B)).LT.1.0D-15) THEN
              M=INT(C-A-B)
              CALL GAMMAX(A,GA)
              CALL GAMMAX(B,GB)
              CALL GAMMAX(C,GC)
              CALL GAMMAX(A+M,GAM)
              CALL GAMMAX(B+M,GBM)
              CALL PSI(A,PA)
              CALL PSI(B,PB)
              IF (M.NE.0) GM=1.0D0
              DO 30 J=1,ABS(M)-1
30               GM=GM*J
              RM=1.0D0
              DO 35 J=1,ABS(M)
35               RM=RM*J
              F0=1.0D0
              R0=1.0D0
              R1=1.0D0
              SP0=0.D0
              SP=0.0D0
              IF (M.GE.0) THEN
                 C0=GM*GC/(GAM*GBM)
                 C1=-GC*(X-1.0D0)**M/(GA*GB*RM)
                 DO 40 K=1,M-1
                    R0=R0*(A+K-1.0D0)*(B+K-1.0)/(K*(K-M))*(1.0-X)
40                  F0=F0+R0
                 DO 45 K=1,M
45                  SP0=SP0+1.0D0/(A+K-1.0)+1.0/(B+K-1.0)-1.0/K
                 F1=PA+PB+SP0+2.0D0*EL+DLOG(1.0D0-X)
                 DO 55 K=1,250
                    SP=SP+(1.0D0-A)/(K*(A+K-1.0))+(1.0-B)/(K*(B+K-1.0))
                    SM=0.0D0
                    DO 50 J=1,M
50                     SM=SM+(1.0D0-A)/((J+K)*(A+J+K-1.0))+1.0/
     &                    (B+J+K-1.0)
                    RP=PA+PB+2.0D0*EL+SP+SM+DLOG(1.0D0-X)
                    R1=R1*(A+M+K-1.0D0)*(B+M+K-1.0)/(K*(M+K))*(1.0-X)
                    F1=F1+R1*RP
                    IF (DABS(F1-HW).LT.DABS(F1)*EPS) GO TO 60
55                  HW=F1
60               HF=F0*C0+F1*C1
              ELSE IF (M.LT.0) THEN
                 M=-M
                 C0=GM*GC/(GA*GB*(1.0D0-X)**M)
                 C1=-(-1)**M*GC/(GAM*GBM*RM)
                 DO 65 K=1,M-1
                    R0=R0*(A-M+K-1.0D0)*(B-M+K-1.0)/(K*(K-M))*(1.0-X)
65                  F0=F0+R0
                 DO 70 K=1,M
70                  SP0=SP0+1.0D0/K
                 F1=PA+PB-SP0+2.0D0*EL+DLOG(1.0D0-X)
                 DO 80 K=1,250
                    SP=SP+(1.0D0-A)/(K*(A+K-1.0))+(1.0-B)/(K*(B+K-1.0))
                    SM=0.0D0
                    DO 75 J=1,M
75                     SM=SM+1.0D0/(J+K)
                    RP=PA+PB+2.0D0*EL+SP-SM+DLOG(1.0D0-X)
                    R1=R1*(A+K-1.0D0)*(B+K-1.0)/(K*(M+K))*(1.0-X)
                    F1=F1+R1*RP
                    IF (DABS(F1-HW).LT.DABS(F1)*EPS) GO TO 85
80                  HW=F1
85               HF=F0*C0+F1*C1
              ENDIF
           ELSE
              CALL GAMMAX(A,GA)
              CALL GAMMAX(B,GB)
              CALL GAMMAX(C,GC)
              CALL GAMMAX(C-A,GCA)
              CALL GAMMAX(C-B,GCB)
              CALL GAMMAX(C-A-B,GCAB)
              CALL GAMMAX(A+B-C,GABC)
              C0=GC*GCAB/(GCA*GCB)
              C1=GC*GABC/(GA*GB)*(1.0D0-X)**(C-A-B)
              HF=0.0D0
              R0=C0
              R1=C1
              DO 90 K=1,250
                 R0=R0*(A+K-1.0D0)*(B+K-1.0)/(K*(A+B-C+K))*(1.0-X)
                 R1=R1*(C-A+K-1.0D0)*(C-B+K-1.0)/(K*(C-A-B+K))
     &              *(1.0-X)
                 HF=HF+R0+R1
                 IF (DABS(HF-HW).LT.DABS(HF)*EPS) GO TO 95
90               HW=HF
95            HF=HF+C0+C1
           ENDIF
        ELSE
           A0=1.0D0
           IF (C.GT.A.AND.C.LT.2.0D0*A.AND.
     &         C.GT.B.AND.C.LT.2.0D0*B) THEN
              A0=(1.0D0-X)**(C-A-B)
              A=C-A
              B=C-B
           ENDIF
           HF=1.0D0
           R=1.0D0
           DO 100 K=1,250
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X
              HF=HF+R
              IF (DABS(HF-HW).LE.DABS(HF)*EPS) GO TO 105
100           HW=HF
105        HF=A0*HF
        ENDIF
        IF (X1.LT.0.0D0) THEN
           X=X1
           C0=1.0D0/(1.0D0-X)**AA
           HF=C0*HF
        ENDIF
        A=AA
        B=BB
        IF (K.GT.120) WRITE(*,115)
115     FORMAT(1X,'Warning! You should check the accuracy')
        RETURN
        END
c
c ******************************************************* c
c
        SUBROUTINE GAMMAX(X,GA)
c
c ******************************************************* c
C
C       ==================================================
C       Purpose: Compute gamma function â(x)
C       Input :  x  --- Argument of â(x)
C                       ( x is not equal to 0,-1,-2,úúú)
C       Output:  GA --- â(x)
C       ==================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,
     &          -0.6558780715202538D0, -0.420026350340952D-1,
     &          0.1665386113822915D0,-.421977345555443D-1,
     &          -.96219715278770D-2, .72189432466630D-2,
     &          -.11651675918591D-2, -.2152416741149D-3,
     &          .1280502823882D-3, -.201348547807D-4,
     &          -.12504934821D-5, .11330272320D-5,
     &          -.2056338417D-6, .61160950D-8,
     &          .50020075D-8, -.11812746D-8,
     &          .1043427D-9, .77823D-11,
     &          -.36968D-11, .51D-12,
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END
c
c ******************************************************* c
c
        SUBROUTINE PSI(X,PS)
c
c ******************************************************* c
C
C       ======================================
C       Purpose: Compute Psi function
C       Input :  x  --- Argument of psi(x)
C       Output:  PS --- psi(x)
C       ======================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        XA=DABS(X)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        S=0.0D0
        IF (X.EQ.INT(X).AND.X.LE.0.0) THEN
           PS=1.0D+300
           RETURN
        ELSE IF (XA.EQ.INT(XA)) THEN
           N=XA
           DO 10 K=1 ,N-1
10            S=S+1.0D0/K
           PS=-EL+S
        ELSE IF (XA+.5.EQ.INT(XA+.5)) THEN
           N=XA-.5
           DO 20 K=1,N
20            S=S+1.0/(2.0D0*K-1.0D0)
           PS=-EL+2.0D0*S-1.386294361119891D0
        ELSE
           IF (XA.LT.10.0) THEN
              N=10-INT(XA)
              DO 30 K=0,N-1
30               S=S+1.0D0/(XA+K)
              XA=XA+N
           ENDIF
           X2=1.0D0/(XA*XA)
           A1=-.8333333333333D-01
           A2=.83333333333333333D-02
           A3=-.39682539682539683D-02
           A4=.41666666666666667D-02
           A5=-.75757575757575758D-02
           A6=.21092796092796093D-01
           A7=-.83333333333333333D-01
           A8=.4432598039215686D0
           PS=DLOG(XA)-.5D0/XA+X2*(((((((A8*X2+A7)*X2+
     &        A6)*X2+A5)*X2+A4)*X2+A3)*X2+A2)*X2+A1)
           PS=PS-S
        ENDIF
        IF (X.LT.0.0) PS=PS-PI*DCOS(PI*X)/DSIN(PI*X)-1.0D0/X
        RETURN
        END
c
c ************************************************************
c
      SUBROUTINE qromo(func,a,b,ss,choose)
c
c ************************************************************
c
      implicit none
      INTEGER JMAX,JMAXP,K,KM
      REAL*8 a,b,func,ss,EPS
      EXTERNAL func,choose
      PARAMETER (EPS=1.d-6, JMAX=14, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint
      INTEGER j
      REAL*8 dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do 11 j=1,JMAX
        call choose(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=h(j)/9.
11    continue
      pause 'too many steps in qromo'
      END
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
C ************************************************************
C
      SUBROUTINE midpnt(func,a,b,s,n)
C
C ************************************************************
C
      implicit none
      INTEGER n
      REAL*8 a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL*8 ddel,del,sum,tnm,x
      if (n.eq.1) then
        s=(b-a)*func(0.5*(a+b))
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      END
C
C ************************************************************
C
      SUBROUTINE midinf(funk,aa,bb,s,n)
C
C ************************************************************
C
      implicit none
      INTEGER n
      REAL*8 aa,bb,s,funk
      EXTERNAL funk
      INTEGER it,j
      REAL*8 a,b,ddel,del,sum,tnm,func,x
      func(x)=funk(1./x)/x**2
      b=1./aa
      a=1./bb
      if (n.eq.1) then
        s=(b-a)*func(0.5*(a+b))
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
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
c ************************************************************
c
      FUNCTION gammp(a,x)
c
c ************************************************************
c
      REAL*8 a,gammp,x
CU    USES gcf,gser
      REAL*8 gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.)pause 'bad arguments in gammp'
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammp=gamser
      else
        call gcf(gammcf,a,x,gln)
        gammp=1.-gammcf
      endif
      return
      END
c
c
c ************************************************************
c
      SUBROUTINE gser(gamser,a,x,gln)
c
c ************************************************************
c
      INTEGER ITMAX
      REAL*8 a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.d-7)
CU    USES gammln
      INTEGER n
      REAL*8 ap,del,sum,gammln
      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.)pause 'x < 0 in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gser'
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      END
c
c ************************************************************

      SUBROUTINE gcf(gammcf,a,x,gln)
c
c ************************************************************
c
      INTEGER ITMAX
      REAL*8 a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.d-7,FPMIN=1.d-30)
CU    USES gammln
      INTEGER i
      REAL*8 an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gcf'
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      END
c
c ************************************************************
c
      FUNCTION gammln(xx)
c
c ************************************************************
c
      REAL*8 gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
c
