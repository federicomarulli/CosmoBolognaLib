c ******************************************************* c
c                                                         c
c  ///     (Standard) Perturbation Theory for         /// c
c  ///    Non-linver Evolution of Power Spectrum      /// c
c  /// in redshift space using Gaussian quadratures   /// c
c                                                         c  
c               Time-stamp: <2011-06-15 19:09:50 ataruya> c   
c                                                         c  
c                                                         c  
c      This code provides monopole (l=0), quadrupole      c  
c      (l=2), and hexadecapole (l=4) moments of           c  
c      redshift-space power spectrum from standard        c  
c      1-loop PT for a given redshift in wCDM cosmology.  c
c      The data format of the output file,                c
c      'pkred_stdPT.dat', is as follows:                  c  
c                                                         c  
c    k, P_lin , P0_EH, P0_PT, P2_EH, P2_PT, P4_EH, P2_PT  c  
c                                                         c  
c       where                                             c
c                                                         c  
c       k:  wavenumber in unit of h/Mpc                   c  
c                                                         c  
c       P_lin: linear matter P(k)  (input power spectrum) c  
c                                                         c  
c       P0_EH, P2_EH, P4_EH : monopole, quadrupole, and   c 
c       hexadecapole moments of linear redshift-space     c
c       P(k) from the no-wiggle transfer function by      c
c       Eisenstein & Hu (1998). Only the linear Kaiser    c
c       effect is taken into account for linear redshift  c  
c       P(k).                                             c  
c                                                         c  
c       P0_PT, P2_PT, P4_PT : monopole, quadrupole, and   c 
c       hexadecapole moments of redshift-space P(k) from  c
c       one-loop standard PT.                             c  
c                                                         c  
c       Note that the units of power spectrum are all     c  
c       Mpc^3/h^3.                                        c  
c                                                         c  
c ******************************************************* c
c
c     Note--. 
c
c     All the parameter 'ikmax' in this program must be 
c     the same. 
c
      program stdPTred
c
      implicit none
c
      integer  ik, ikmax, ik_max, ik_max1, ik_max2
      integer  itype, ibox, iparams
      parameter(ikmax=3000)
      real*8  ak(ikmax), pk(ikmax)
      real*8  akPT(ikmax), pk0PT(ikmax), pk2PT(ikmax), pk4PT(ikmax)
      real*8  akEH(ikmax), pk0EH(ikmax), pk2EH(ikmax), pk4EH(ikmax)
      real*8  h, Tcmb, n_s, sigma8, Omega_m, Omega_b, Omega_v, w_de
      real*8   pi, Lbox, zred
      character infile*50
      common /cosm_param/ Omega_b, Omega_m, h, Tcmb, n_s, w_de
      common /pk_data/ ak, pk, ik_max
      common /pkPT/ akPT, pk0PT, pk2PT, pk4PT, ik_max1
      common /pkEH/ akEH, pk0EH, pk2EH, pk4EH, ik_max2
      pi = 4.d0 * atan(1.d0)
c     -------------------------------------------------
c
c     Cosmological parameters for linear growh factor (default values)
c
c     Nishimichi et al. (2009)
c
      h = 0.701d0
      Tcmb = 2.726d0
      n_s = 0.96d0
      sigma8 = 0.817d0
      Omega_m = 0.279d0
      Omega_b = 0.165 * Omega_m
      Omega_v = 1.d0-Omega_m 
      w_de = -1.0
c
      write(6,*)
      write(6,*) '*** PT calculation of P(k) in redshift space ***'
      write(6,*)
      write(6,*)
      write(6,*) 'Output redshift ?'
      read(5,*) zred
      write(6,*) 'including finite-volume effect? :type (0 or 1, L_box)'
      write(6,*) ' (e.g., (0, 300) for boxsize L_box=300Mpc/h ) '
      write(6,*) ' (      (1, ???) for ignoring finite-volume effect) '
      read(5,*) ibox, Lbox
c
c ///// Set cosmological parameters ///// c
c
 5    write(6,*) '[1] omega_m =',Omega_m
      write(6,*) '[2] omega_b =',Omega_b
      write(6,*) '[3]       h =',h
      write(6,*) '[4]     n_s =',n_s
      write(6,*) '[5]  sigma8 =',sigma8
      write(6,*) '[6]       w =',w_de
      write(6,*) 
      write(6,*) 'Note--. spatial curvature is forced to be flat,'
      write(6,*) '        i.e., omega_m = 1 - omega_DE           '
c
      write(6,*)
      write(6,*) 'change cosmological parameter? [1-6] or n[0]'
      read(5,*)  iparams
      if(iparams.eq.0) goto 8
      if(iparams.eq.1) then
         write(6,*) 'type omega_m'
         read(5,*) Omega_m
         Omega_v = 1.d0 - Omega_m
      elseif(iparams.eq.2) then
         write(6,*) 'type omega_b'
         read(5,*) Omega_b
      elseif(iparams.eq.3) then
         write(6,*) 'type h'
         read(5,*) h
      elseif(iparams.eq.4) then
         write(6,*) 'type n_s'
         read(5,*) n_s
      elseif(iparams.eq.5) then
         write(6,*) 'type sigma8'
         read(5,*) sigma8
      elseif(iparams.eq.6) then
         write(6,*) 'type w'
         read(5,*) w_de
      endif
      goto 5
c
c     /////// Load transfer function ///////
c
 8    call load_transfer_function
c
      write(6,*) ' loading (linear) matter power spectrum, done '
c
c     /////// Sigma8 normalization ///////
c
      call normalization_trapez(sigma8)
c
      write(6,*) ' normalization by sigma8 done '
c
c     //////// Truncation of low-k, high-k  in linear P(k) ////////
c
      call truncation_k(ibox, Lbox)
c
c     /////// 1-loop correction of P(k) ///////
c
      call calc_one_loop_pk
c
      write(6,*) ' 1-loop P(k) done '
c
c     /////// Summing up all contributions ///////
c
      call calc_pkred(zred,sigma8)
c
      write(6,*) ' summing up all the contributions done '
c
c     /////// Save output data ///////
c
      open(10,file='pkred_stdPT.dat',status='unknown')
c
      do ik =1, ik_max
         write(10,'(1p8e18.10)') ak(ik), pk(ik), 
     &        pk0EH(ik), pk0PT(ik), pk2EH(ik), 
     &        pk2PT(ik), pk4EH(ik), pk4PT(ik) 
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
      common /pk_data/ ak, pk, ik_max
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
      common /pk_data/ ak, pk, ik_max
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
      kmin = 5.d-4   ! default value
c
      if(ibox.eq.0) kmin = 2.d0 * pi / Lbox
c
      kmax = 10.d0    ! default value
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
      subroutine calc_one_loop_pk
c
c ******************************************************* c
c
      implicit none
      integer ik, ik_max, ikmax, isub
      integer ix, ixmax13, ixmax22
      parameter(ikmax=3000)
      parameter(ixmax13=100, ixmax22=600)
      real*8  pi
      real*8  ak(ikmax), pk(ikmax)
      real*8  pk13_B00(ikmax), pk13_B11(ikmax), pk13_B12(ikmax)
      real*8  pk13_B22(ikmax), pk13_B23(ikmax), pk22_A00(ikmax)
      real*8  pk22_A11(ikmax), pk22_A12(ikmax), pk22_A22(ikmax)
      real*8  pk22_A23(ikmax), pk22_A24(ikmax), pk22_A33(ikmax)
      real*8  pk22_A34(ikmax), pk22_A44(ikmax)
      real*8  kmin, kmax, xmin, xmax, mumin, mumax
      real*8  k, pklink
      real*8  w13(ixmax13), x13(ixmax13), w22(ixmax22), x22(ixmax22)
      real*8  fp13, integ_fp22
      common /pk_data/ ak, pk, ik_max
      common /wave_number/  k, xmin, xmax
      common /PT_pk_B/ pk13_B00, pk13_B11, pk13_B12, pk13_B22, pk13_B23
      common /PT_pk_A/ pk22_A00, pk22_A11, pk22_A12, pk22_A22, pk22_A23,
     &     pk22_A24, pk22_A33, pk22_A34, pk22_A44      
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
         pk13_B00(ik)= 0.d0
         pk13_B11(ik)= 0.d0
         pk13_B12(ik)= 0.d0
         pk13_B22(ik)= 0.d0
         pk13_B23(ik)= 0.d0
         pk22_A00(ik)= 0.d0
         pk22_A11(ik)= 0.d0
         pk22_A12(ik)= 0.d0
         pk22_A22(ik)= 0.d0
         pk22_A23(ik)= 0.d0
         pk22_A24(ik)= 0.d0
         pk22_A33(ik)= 0.d0
         pk22_A34(ik)= 0.d0
         pk22_A44(ik)= 0.d0
c
         xmin = kmin / k
         xmax = kmax / k
c
c     ////// Gauss-Legendre integration of p13(k), p22(k) //////  c
c
         call gauleg(log(xmin),log(xmax),x13,w13,ixmax13)
c
         do ix=1, ixmax13
            x13(ix)= dexp(x13(ix))
            pk13_B00(ik) = pk13_B00(ik) + w13(ix) * fp13(1, x13(ix))
            pk13_B12(ik) = pk13_B12(ik) + w13(ix) * fp13(2, x13(ix))
            pk13_B22(ik) = pk13_B22(ik) + w13(ix) * fp13(3, x13(ix))
            pk13_B23(ik) = pk13_B23(ik) + w13(ix) * fp13(4, x13(ix))
         enddo
c
         if(k.lt.0.2) isub =200 
         if(k.ge.0.2) isub =0 
c
            call gauleg(log(xmin),log(xmax),x22,w22,ixmax22-isub)
c
            do ix=1, ixmax22-isub
               x22(ix)= dexp(x22(ix))
               pk22_A00(ik) = pk22_A00(ik)+w22(ix)*integ_fp22(1,x22(ix))
               pk22_A12(ik) = pk22_A12(ik)+w22(ix)*integ_fp22(2,x22(ix))
               pk22_A22(ik) = pk22_A22(ik)+w22(ix)*integ_fp22(3,x22(ix))
               pk22_A24(ik) = pk22_A24(ik)+w22(ix)*integ_fp22(4,x22(ix))
               pk22_A33(ik) = pk22_A33(ik)+w22(ix)*integ_fp22(5,x22(ix))
               pk22_A34(ik) = pk22_A34(ik)+w22(ix)*integ_fp22(6,x22(ix))
               pk22_A44(ik) = pk22_A44(ik)+w22(ix)*integ_fp22(7,x22(ix))
            enddo
c
         pk13_B00(ik) = pk13_B00(ik) * pklink * k**3 / (2.*pi)**2
         pk13_B11(ik) = pk13_B00(ik) * 3.d0 
         pk13_B12(ik) = pk13_B12(ik) * pklink * k**3 / (2.*pi)**2
         pk13_B22(ik) = pk13_B22(ik) * pklink * k**3 / (2.*pi)**2
         pk13_B23(ik) = pk13_B23(ik) * pklink * k**3 / (2.*pi)**2
         pk22_A00(ik) = 2.d0 * pk22_A00(ik) * k**3 / (2.*pi)**2
         pk22_A11(ik) = pk22_A00(ik) * 4.d0
         pk22_A12(ik) = 2.d0 * pk22_A12(ik) * k**3 / (2.*pi)**2
         pk22_A22(ik) = 2.d0 * pk22_A22(ik) * k**3 / (2.*pi)**2
         pk22_A23(ik) = pk22_A12(ik) * 2.d0
         pk22_A24(ik) = 2.d0 * pk22_A24(ik) * k**3 / (2.*pi)**2
         pk22_A33(ik) = 2.d0 * pk22_A33(ik) * k**3 / (2.*pi)**2
         pk22_A34(ik) = 2.d0 * pk22_A34(ik) * k**3 / (2.*pi)**2
         pk22_A44(ik) = 2.d0 * pk22_A44(ik) * k**3 / (2.*pi)**2
c
      write(6,'(i3,1p4e18.10)') ik,k,pklink,pk13_B00(ik),pk22_A00(ik)
c
 10   continue
c
      end
c
c ******************************************************* c
c
      function  fp13(ip, x)
c
c ******************************************************* c
c
c     ip=1 : kernel for pk13_B00
c     ip=2 : kernel for pk13_B12
c     ip=3 : kernel for pk13_B22
c     ip=4 : kernel for pk13_B23
c
      implicit none
      integer ip
      real*8 fp13, x, k, s, pklin
      real*8 xmin, xmax
      common /wave_number/  k, xmin, xmax
c
      fp13 = 0.d0
c
      if(x.lt.0.) then
         write(6,*) ' warning: fp13 negative !! : x=',x
         fp13 = 0.d0
      elseif(x .lt. 5.d-3) then
         if(ip.eq.1) fp13 = -2.d0/3.d0 + 232.d0/315.d0*x*x 
     &        - 376.d0/735.d0*x*x*x*x
         if(ip.eq.2) fp13 = -2.d0/3.d0 - 32.d0/35.d0*x*x 
     &        + 96.d0/245.d0*x*x*x*x
         if(ip.eq.3) fp13 = -4.d0/3.d0 + 48.d0/35.d0*x*x 
     &        - 48.d0/49.d0*x*x*x*x
         if(ip.eq.4) fp13 = -2.d0/3.d0
      elseif(x .gt. 5.d2) then
         s = 1.d0 / x
         if(ip.eq.1) fp13 = -122.d0/315.d0 + 8.d0/105.d0*s*s
     &        -40.d0/1323.d0*s*s*s*s
         if(ip.eq.2) fp13 = -166.d0/105.d0 + 96.d0/245.d0*s*s
     &        -32.d0/735.d0*s*s*s*s
         if(ip.eq.3) fp13 = -92.d0/105.d0 + 48.d0/245.d0*s*s
     &        -16.d0/245.d0*s*s*s*s
         if(ip.eq.4) fp13 = -2.d0/3.d0
      elseif(x.ge.0.995 .and. x.le.1.005) then
         if(ip.eq.1) fp13 = (-22.d0 + 2.d0*(x-1.d0) - 29.d0*(x-1.d0)**2)
     &        / 63.d0
         if(ip.eq.2) fp13 = (-26.d0 - 12.d0*( (x-1.d0) - (x-1.d0)**2 ) )
     &        / 21.d0
         if(ip.eq.3) fp13 = (-16.d0 - 18.d0* (x-1.d0)**2 ) / 21.d0
         if(ip.eq.4) fp13 = -2.d0/3.d0
      else 
         if(ip.eq.1) fp13 = ( 12./x/x - 158. + 100.*x*x -42.*x*x*x*x + 
     &        3./(x*x*x) * (x*x-1.)**3 * (7.*x*x + 2.) * 
     &        dlog(abs((x+1.)/(x-1.))) )/252.
         if(ip.eq.2) fp13 = ( 18./x/x - 178. - 66.*x*x +18.*x*x*x*x - 
     &        9./(x*x*x) * (x*x-1.)**4 * 
     &        dlog(abs((x+1.)/(x-1.))) )/168.d0
         if(ip.eq.3) fp13 = ( 18./x/x - 218. + 126.*x*x -54.*x*x*x*x + 
     &        9./(x*x*x) * (x*x-1.)**3 * (3.*x*x+1.) * 
     &        dlog(abs((x+1.)/(x-1.))) )/168.d0
         if(ip.eq.4) fp13 = -2.d0/3.d0
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
      function integ_fp22(ip, x)
c
c ******************************************************* c
c
      implicit none
      integer ip, imu, imu_max
      parameter(imu_max=10)
      real*8  integ_fp22, xmin, xmax, mumin, mumax, fp22
      real*8  k, x, mu, wmu(imu_max), mmu(imu_max)
      common /wave_number/  k, xmin, xmax
c     ------------------------------------------  c
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
         integ_fp22 = integ_fp22 + wmu(imu) * fp22(ip, x, mmu(imu))
      enddo
c
      end
c
c ******************************************************* c
c
      function  fp22(ip, x, mu)
c
c ******************************************************* c
c
c     ip=1 for kernel of pk22_A00 
c     ip=2 for kernel of pk22_A12
c     ip=3 for kernel of pk22_A22
c     ip=4 for kernel of pk22_A24
c     ip=5 for kernel of pk22_A33
c     ip=6 for kernel of pk22_A34
c     ip=7 for kernel of pk22_A44
c
      implicit none
      integer ip
      real*8 fp22
      real*8 mu, k, x, pklin1, pklin2
      real*8 mumin, mumax, xmax, xmin
      common /wave_number/  k, xmin, xmax
c     --------------------------------------- c
c
      if(ip.eq.1) then
         fp22 = (3.*x + 7.*mu - 10.*mu*mu*x)**2 / 98.
      elseif(ip.eq.2) then
         fp22 = (1.d0-mu*mu) * (7.-6.*x*x-42.*mu*x+48.*(x*mu)**2) / 28.
      elseif(ip.eq.3) then
         fp22 = ( -49. + 637.*mu*mu + 42.*x*mu*( 17. - 45.*mu*mu ) 
     &        + 6.*x*x*( 19. - 157.*mu*mu + 236.*mu**4) ) / 196.
      elseif(ip.eq.4) then
         fp22 = 3. * ( x*(1.-mu*mu) )**2 / 16.
      elseif(ip.eq.5) then
         fp22 = ( -7. + 35.*mu*mu + 54.*x*mu - 110.*x*mu*mu*mu
     &        + 6.*x*x - 66.*(x*mu)**2 + 88.*(x*mu*mu)**2 ) / 14.
      elseif(ip.eq.6) then
         fp22 = (1.d0-mu*mu) * (2.-3.*x*x-12.*mu*x+15.*(x*mu)**2) / 8.
      elseif(ip.eq.7) then
         fp22 = ( -4. + 12.*mu*mu + 24.*x*mu - 40.*x*mu*mu*mu
     &        + 3.*x*x - 30.*(x*mu)**2 + 35.*(x*mu*mu)**2 ) / 16.
      endif
c
      call find_pk(k*x, pklin1)
      call find_pk(k*sqrt(1.+x*x-2.*mu*x), pklin2)
c
      fp22 = fp22 * x * pklin1 * pklin2 / (1.+x*x-2.*mu*x)**2
c
      end
c
c ******************************************************* c
c
      subroutine calc_pkred(zred, sigma8)
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
      real*8 zred, sigma8, growth, f, dd2, ff
      real*8  Omega_b, Omega_m, h, Tcmb, n_s, w_de
      real*8  ak(ikmax), pk(ikmax), pk_EH(ikmax)
      real*8  pk13_B00(ikmax), pk13_B11(ikmax), pk13_B12(ikmax)
      real*8  pk13_B22(ikmax), pk13_B23(ikmax), pk22_A00(ikmax)
      real*8  pk22_A11(ikmax), pk22_A12(ikmax), pk22_A22(ikmax)
      real*8  pk22_A23(ikmax), pk22_A24(ikmax), pk22_A33(ikmax)
      real*8  pk22_A34(ikmax), pk22_A44(ikmax)
      real*8  akPT(ikmax), pk0PT(ikmax), pk2PT(ikmax), pk4PT(ikmax)
      real*8  akEH(ikmax), pk0EH(ikmax), pk2EH(ikmax), pk4EH(ikmax)
      real*8  factA, factB
      common /growth_index/ ff 
      common /cosm_param/ Omega_b, Omega_m, h, Tcmb, n_s, w_de
      common /pk_data/ ak, pk, ik_max
      common /PT_pk_B/ pk13_B00, pk13_B11, pk13_B12, pk13_B22, pk13_B23
      common /PT_pk_A/ pk22_A00, pk22_A11, pk22_A12, pk22_A22, pk22_A23,
     &     pk22_A24, pk22_A33, pk22_A34, pk22_A44
      common /pkPT/ akPT, pk0PT, pk2PT, pk4PT, ik_max1
      common /pkEH/ akEH, pk0EH, pk2EH, pk4EH, ik_max2
c     -------------------------------------------------------- c
      ik_max1 = ik_max
      ik_max2 = ik_max
c
c     growth factor and its logarithmic derivative
c     (assuming flat universe)
      dd2 = growth(zred, Omega_m, 1.d0-Omega_m, w_de)**2
      ff = f(zred, Omega_m, 1.d0-Omega_m, w_de)
c
c     ////// one-loop power spectrum  ////// c
c
      do ik=1, ik_max1
         akPT(ik) = ak(ik)
         pk0PT(ik) = (1.d0+2.d0/3.d0*ff+1.d0/5.d0*ff*ff) * pk(ik)
         pk2PT(ik) = ( 4.d0/3.d0*ff + 4.d0/7.d0*ff*ff ) * pk(ik) 
         pk4PT(ik) = 8.d0/35.d0*ff*ff * pk(ik) 
c
         pk0PT(ik) = pk0PT(ik) + dd2 * ( 
     &        factA(0,0) * pk22_A00(ik) + 
     &        factA(1,0) * ( ff*pk22_A11(ik) + ff**2*pk22_A12(ik) ) + 
     &        factA(2,0) * ( ff**2*pk22_A22(ik) + ff**3*pk22_A23(ik) +
     &        ff**4*pk22_A24(ik) ) + 
     &        factA(3,0) * ( ff**3*pk22_A33(ik) + ff**4*pk22_A34(ik) ) +
     &        factA(4,0) * ff**4*pk22_A44(ik)  ) 
     &        + dd2 * (
     &        factB(0,0) * pk13_B00(ik) +
     &        factB(1,0) * ( ff*pk13_B11(ik) + ff**2*pk13_B12(ik) ) +
     &        factB(2,0) * ( ff**2*pk13_B22(ik) + ff**3*pk13_B23(ik) ) ) 
c
         pk2PT(ik) = pk2PT(ik) + dd2 * ( 
     &        factA(0,2) * pk22_A00(ik) + 
     &        factA(1,2) * ( ff*pk22_A11(ik) + ff**2*pk22_A12(ik) ) + 
     &        factA(2,2) * ( ff**2*pk22_A22(ik) + ff**3*pk22_A23(ik) +
     &        ff**4*pk22_A24(ik) ) + 
     &        factA(3,2) * ( ff**3*pk22_A33(ik) + ff**4*pk22_A34(ik) ) +
     &        factA(4,2) * ff**4*pk22_A44(ik)  ) 
     &        + dd2 * (
     &        factB(0,2) * pk13_B00(ik) +
     &        factB(1,2) * ( ff*pk13_B11(ik) + ff**2*pk13_B12(ik) ) +
     &        factB(2,2) * ( ff**2*pk13_B22(ik) + ff**3*pk13_B23(ik) ) ) 
c
         pk4PT(ik) = pk4PT(ik) + dd2 * ( 
     &        factA(0,4) * pk22_A00(ik) + 
     &        factA(1,4) * ( ff*pk22_A11(ik) + ff**2*pk22_A12(ik) ) + 
     &        factA(2,4) * ( ff**2*pk22_A22(ik) + ff**3*pk22_A23(ik) +
     &        ff**4*pk22_A24(ik) ) + 
     &        factA(3,4) * ( ff**3*pk22_A33(ik) + ff**4*pk22_A34(ik) ) +
     &        factA(4,4) * ff**4*pk22_A44(ik)  ) 
     &        + dd2 * (
     &        factB(0,4) * pk13_B00(ik) +
     &        factB(1,4) * ( ff*pk13_B11(ik) + ff**2*pk13_B12(ik) ) +
     &        factB(2,4) * ( ff**2*pk13_B22(ik) + ff**3*pk13_B23(ik) ) ) 
c
         pk0PT(ik) = pk0PT(ik) * dd2
         pk2PT(ik) = pk2PT(ik) * dd2
         pk4PT(ik) = pk4PT(ik) * dd2
c
      enddo
c
c     ////// no-wiggle linear power spectrum  ////// c
c
      call no_wiggle_Pk(sigma8, ik_max, ak, pk_EH)
c
      do ik=1, ik_max2
         akEH(ik) = ak(ik)
         pk0EH(ik) = (1.d0+2.d0/3.d0*ff+1.d0/5.d0*ff*ff)* pk_EH(ik)*dd2
         pk2EH(ik) = ( 4.d0/3.d0*ff + 4.d0/7.d0*ff*ff )* pk_EH(ik)*dd2
         pk4EH(ik) = 8.d0/35.d0*ff*ff * pk_EH(ik)*dd2
      enddo
c
      end
c
c ************************************************ c
c
      function factA(n, l)
c
c ************************************************ c
c
c     (2l+1)/2 * integ dmu  mu^(2n) * P_l(mu)
c
      implicit none
      integer n, l
      real*8 factA, nn
c     ---------------------------- c
      nn = dble(n)
c
      if(l.eq.0) factA = 1.d0 / (4.d0*nn + 2.d0)
      if(l.eq.2) factA = 5.d0*nn / (4.d0*nn**2 + 8.d0*nn + 3.d0)
      if(l.eq.4) factA = 18.d0*nn * (nn-1.d0) 
     &     / (8.d0*nn**3 + 36.d0*nn**2 + 46.d0*nn + 15.d0)
c
      factA = factA * (1.d0 + (-1.d0)**(2.*nn)) 
c
      end
c
c ************************************************ c
c
      function factB(n, l)
c
c ************************************************ c
c
c     (2l+1)/2 * integ dmu  (1+f*mu^2) * mu^(2n) * P_l(mu)
c
      implicit none
      integer n, l
      real*8 factB, nn, ff
      common /growth_index/ ff 
c     ---------------------------- c
      nn = dble(n)
c
      if(l.eq.0) then
         factB = 2.d0*nn*(1.d0+ff) + 3.d0 + ff
         factB = factB / (8.d0*nn**2 + 16.d0*nn + 6.d0)
      elseif(l.eq.2) then
         factB = 5.d0 * (2.d0*nn**2*(1.d0+ff) + nn*(5.d0+3.d0*ff) + ff)
         factB = factB  / (1.d0+2.d0*nn) /(3.d0+2.d0*nn)/(5.d0+2.d0*nn)
      elseif(l.eq.4) then
         factB = 2.d0*nn**2*(1.d0+ff) + nn*(5.d0+3.d0*ff) + ff -7.d0 
         factB = factB * 18.d0 * nn / (1.d0+2.d0*nn) / (3.d0+2.d0*nn) 
         factB = factB / (5.d0+2.d0*nn) / (7.d0+2.d0*nn)
      endif
c
      factB = factB * (1.d0 + (-1.d0)**(2.*nn)) 
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
      real*8 omegab, omega0, h, Tcmb, n_s, w_de
      common /cosm_param/ omegab, omega0, h, Tcmb, n_s, w_de
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
      function growth(zred, Omega_m, Omega_v, w_de)
c
c ******************************************************* c
c
c     Linear growth factor for flat cosmology
c
      implicit none
      real*8 growth, zred, Omega_m, Omega_v, w_de, a, b, c
      real*8 zred1, zi, zz, gz, g0
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
      function f(zred, Omega_m, Omega_v, w_de)
c
c ******************************************************* c
c
c     d lnD_+(z) / d ln a,   as function of redshift 
c
      implicit none
c
      real*8  f, zred, zred1, Omega_m, Omega_v, w_de
      real*8  a, b, c, zi, zz, g1, g2
c     ---------------------------------------------------
c
      zred1 = 1.d0 + zred
      zi = - Omega_v / Omega_m
      zz = zi * zred1**(3*w_de) 
c
      a = 1.d0 - 1.d0 / (3.d0 * w_de)
      b = 1.5d0 - 1.d0 / (2.d0 * w_de) 
      c = 2.d0 - 5.d0 / (6.d0 * w_de)
      call HYGFX(a,b,c,zz, g1)  
c
      a = - 1.d0 / (3.d0 * w_de)
      b = 0.5d0 - 1.d0 / (2.d0 * w_de)
      c = 1 - 5.d0/ (6.d0 * w_de)
      call HYGFX(a,b,c,zz, g2)  
c
      f = 1.d0 + 3.d0*(w_de-1.d0)/(6.d0*w_de-5.d0) * zz * g1 / g2
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
