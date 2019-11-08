c ********************************************************** c
c
      program fit_pkred
c
c                  Time-stamp: <2011-06-15 19:12:22 ataruya> c
c ********************************************************** c
c
c     reading the data of non-linear P(k) created from read_pk2.f 
c     and computing redshift P(k) by fitting sigmav2 to N-body data 
c
c     This file uses following data:
c     
c     pkrenorm2_Closure.dat (or pkrenorm2_RPT.dat)
c     pknbody.dat
c
      implicit none
      integer ik, ikmax, ik_max1, ik_max2, ik_max_nbody, ichoice
      integer is, is_max, iparams, iter, itst, j, iguess, iks_max
      integer imodel, imodel_s
      parameter(ikmax=5000)
      real*8  Tcmb, sigma8, n_s, Omega_m, Omega_b, Omega_v, h
      real*8  zred, pk_lin_red
      real*8  ak1(ikmax), ak2(ikmax), pk_lin(ikmax), sigmav2
      real*8  pk_dd(ikmax), pk_dt(ikmax), pk_tt(ikmax)
      real*8  dummy, pk_sim0(ikmax), pk_sim2(ikmax)
      real*8  pk_EH(ikmax), pk_EH0, pk_EH2 , pk_EH4
      real*8  pk_lin0, pk_lin2, pk_lin4, ss0, ss2, ss4
      real*8  aks(ikmax), pk0(ikmax), pk2(ikmax), pk4(ikmax)
      real*8  var0(ikmax), var2(ikmax), kmin, kmax, sigmav2_guess
      real*8  aksim(ikmax), pksim(ikmax), varsim(ikmax)
      real*8  as(ikmax), xi0(ikmax), xi2(ikmax), xi4(ikmax)
c
      integer ia(1)
      external calc_pkred
      real*8  alambda, chisq, a(1), covar(1,1), alpha(1,1), ochisq
      common  /cosm_param/ Omega_b, Omega_m, h, Tcmb, n_s
      common  /pk_data/ ak1, pk_dd, pk_dt, pk_tt, ik_max1
      common  /pk_nbody/ ak2, pk_sim0, pk_sim2, var0, var2, ik_max2
      common  /redshift/ zred
      common  /model_s/ imodel_s
      common /kmin_kmax/  kmin, kmax
      common /pkred_data/ aks, pk0, pk2, pk4, iks_max
c     -------------------------------------------------    c
c
c     Cosmological parameters for linear growh factor (default values)
c
c     Nishimichi
c
      h = 0.701d0
      Tcmb = 2.726d0
      n_s = 0.96d0
      sigma8 = 0.817d0
      Omega_m = 0.279d0
      Omega_b = 0.165 * Omega_m
      Omega_v = 1.d0-Omega_m 
c
c ///// Set cosmological parameters ///// c
c
c
      write(6,*) ' type redshift ' 
      read(5,*) zred
c
 5    write(6,*) '[1] omega_m =',Omega_m
      write(6,*) '[2] omega_b =',Omega_b
      write(6,*) '[3]       h =',h
      write(6,*) '[4]     n_s =',n_s
      write(6,*) '[5]  sigma8 =',sigma8
      write(6,*) 
      write(6,*) 'Note--. spatial curvature is forced to be flat'
      write(6,*)
      write(6,*) 'change cosmological parameter? [1-5] or n[0]'
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
      endif
      goto 5
c
c     ////////// Read data files //////////
c
 8    write(6,*) ' fit N-body to SPT[0] or CLA[1] ? '
      read(5,*) imodel
      if(imodel.eq.0) 
     &     open(9,file='pkstd2.dat',status='unknown') 
      if(imodel.eq.1) 
     &     open(9,file='pkrenorm2_Closure.dat',status='unknown') 
      open(10,file='pknbody_red.dat',status='unknown') 
c
      do ik = 1, ikmax
      read(9,*,END=10) ak1(ik), pk_EH(ik), pk_lin(ik), pk_dd(ik),
     &        pk_dt(ik), pk_tt(ik)
      enddo
c
 10   ik_max1 = ik - 1
c
      write(6,'(a,1p1e12.6,a,1p1e12.6)') 
     &     'kmin_model=',ak1(1),', kmax_model=',ak1(ik_max1)
c
      do ik = 1, ikmax
      read(10,*,END=20) ak2(ik), dummy, pk_sim0(ik), var0(ik),
     &        dummy, pk_sim2(ik), var2(ik), dummy, dummy, dummy
ccccc         read(10,*,END=20) ak2(ik), pk_dummy, pk_sim0(ik), var0(ik)
      enddo
c
 20   ik_max2 = ik - 1 
c
      write(6,'(a,1p1e12.6,a,1p1e12.6)') 
     &     'kmin_sim=',ak2(1),', kmax_sim=',ak2(ik_max2)
c
      close(9)
      close(10)
c
c     ////////// Fitting sigmav2 //////////
c
c     note--. The actual fitting parameter is sigmav, not sigmav2 
c
      write(6,*)
      write(6,*) 'type fitting range kmin, kmax '
      read(5,*) kmin, kmax
      if(kmin.lt.ak2(1) .or. kmax.gt.ak2(ik_max2-1)) then
         write(6,*) ' out of range !! '
         stop
      endif
c
      write(6,*)
      write(6,*) ' Select model of redshift-space distortion: '
      write(6,*) '         Scoccimarro (2004)              [0]'
      write(6,*) '         Kaiser + Gaussian-FoG           [1]'
      write(6,*) '         Scoccimarro (2004) with exp-FoG [2]'
      read(5,*)  imodel_s
      if(imodel_s.lt.0 .or. imodel_s.gt.2) stop 
c
      write(6,*)
      write(6,*) 'choose data set: pk0 only [0], pk0 & pk2 [1]' 
      read(5,*) ichoice
c
      call resorting_nbody_data(ichoice, kmin, kmax, aksim, pksim, 
     &     varsim, ik_max_nbody)
c
      call get_initial_sigmav2(sigmav2_guess)
c
      write(6,*) 'sigmav2=', sigmav2_guess
c
      write(6,*) 'use this value for initial guess ? y[0], n[1]'
      read(5,*) iguess
      if(iguess.eq.1) then
         write(6,*) ' input sigmav2(initial) '
         read(5,*)  sigmav2_guess
      endif
c
      ia(1) = 1
      a(1) = dsqrt(sigmav2_guess)
c
      alambda = -1.d0
      call mrqmin(aksim, pksim, varsim, ik_max_nbody, a, ia, 1, 
     &     covar, alpha, 1, chisq, calc_pkred, alambda)
c
      j = 1
      itst = 0
c
 99   write(6,'(/1x,a,i2,t18,a,g10.4,t43,a,g9.2)') 'Iteration #',j,
     *     'Chi-squared:',chisq,'ALAMDA:',alambda
      write(6,'(A,g12.4)') 'sigmav2=',a(1)**2
      j = j + 1
      ochisq = chisq
      call mrqmin(aksim, pksim, varsim, ik_max_nbody, a, ia, 1, 
     &     covar, alpha, 1, chisq, calc_pkred, alambda)
      if(chisq.gt.ochisq) then 
         itst=0
      elseif (abs(ochisq-chisq).lt.0.1) then
         itst=itst+1
      endif
      if (itst.lt.10) then
         goto 99
      endif
      alambda=0.d0
      call mrqmin(aksim, pksim, varsim, ik_max_nbody, a, ia, 1, 
     &     covar, alpha, 1, chisq, calc_pkred, alambda)
c
ccc         write(*,*) 'Uncertainties:'
ccc         write(*,*) sqrt(covar(1,1))
c
      write(6,*)
      write(6,*) 'Final result:    '
      write(6,'(A,1p1g16.6,A,1p1g16.6,A,1p1g16.6)')  
     &     'sigmav2=',a(1)**2, ', chi^2=',chisq, 
     &     ', d_sigmav', sqrt(covar(1,1))
      write(6,*)
      write(6,*) '(c.f.) sigmav2_init=', sigmav2_guess
c
c
c     ////////// Output redshift P(k) data //////////
c
      if(imodel.eq.0)
     &     open(11,file='pkstd_red_sigmav2_fit.dat',
     &     status='unknown')
      if(imodel.eq.1)
     &     open(11,file='pkrenorm_Closure_red_sigmav2_fit.dat',
     &     status='unknown')
c
      sigmav2 = a(1)**2
c
      iks_max = ik_max1
      do ik = 1, ik_max1
c
         pk_EH0 = pk_lin_red(0, pk_EH(ik))
         pk_EH2 = pk_lin_red(2, pk_EH(ik))
         pk_EH4 = pk_lin_red(4, pk_EH(ik))
         pk_lin0 = pk_lin_red(0, pk_lin(ik))
         pk_lin2 = pk_lin_red(2, pk_lin(ik))
         pk_lin4 = pk_lin_red(4, pk_lin(ik))
c
         aks(ik) = ak1(ik)
         call pk_red(ik, sigmav2, ss0, ss2, ss4)
         pk0(ik) = ss0
         pk2(ik) = ss2
         pk4(ik) = ss4
c
         write(11,'(1p10e18.10)') aks(ik), pk_EH0, 
     &        pk_lin0, pk0(ik), pk_EH2, pk_lin2, pk2(ik), 
     &        pk_EH4, pk_lin4, pk4(ik)
c         
      enddo
c
      write(6,*)
      write(6,*) ' Save fitting data to   '
      if(imodel.eq.0)
     &     write(6,*) '     pkstd_red_sigmav2_fit.dat '
      if(imodel.eq.1)
     &     write(6,*) '     pkrenorm_Closure_red_sigmav2_fit.dat '
c
      close(11)
c
c     ////////// Output redshift xi(s) data //////////
c
      if(imodel.eq.1) then
         write(6,*) ' Save Xi_0(s), Xi_2(s) and Xi_4(s) ? Yes[0], No[1]'
         read(5,*) iparams
         if(iparams.eq.1) then
            write(6,*) ' Bye ! '
            stop
         else
            open(10,file='xi_red_sigmav2_fit.dat',status='unknown')
c     
            kmin = aks(1)
            kmax = aks(iks_max)
            call calc_xi(as, xi0, xi2, xi4, is_max)
c
            do is = 1, is_max
               write(10,'(1p4e18.10)') as(is), xi0(is), xi2(is), xi4(is)
            enddo
            close(10)
            write(6,*) ' Save xi_red_sigmav2_fit.dat'
c
         endif
      endif
c
      end
c
c ******************************************************* c
c
      subroutine calc_xi(as, xi0, xi2, xi4, is_max)
c
c ******************************************************* c
c
      implicit none
      integer  ikmax, is, ismax, is_max, ix, ixmax
      parameter(ikmax=5000, ismax=501, ixmax=1000)
      real*8  as(ikmax), xi0(ikmax), xi2(ikmax), xi4(ikmax)
      real*8  smin, smax, pi, j0, j2, j4, kmin, kmax
      real*8  pk0, pk2, pk4
      parameter(smin=5.d0, smax=200.d0) 
      real*8  k, kk(ixmax), ww(ixmax)
      common /kmin_kmax/  kmin, kmax
      pi = 4.d0 * datan(1.d0)
c     -------------------------------------------------
c
      is_max = ismax
c
      call gauleg(dlog(kmin), dlog(kmax), kk, ww, ixmax)
c
      do is = 1, ismax
c
         as(is) =  smin * (smax/smin) ** (dble(is-1)/dble(ismax-1))
c
         xi0(is) = 0.d0
         xi2(is) = 0.d0
         xi4(is) = 0.d0
c
         do ix=1, ixmax
c
            k = dexp(kk(ix)) 
            call find_pkred(1, k, pk0)
            call find_pkred(2, k, pk2)
            call find_pkred(3, k, pk4)
c
            xi0(is) = xi0(is) + ww(ix) * pk0 * j0(as(is)*k) 
     &           * k**3 / (2.d0*pi*pi)
            xi2(is) = xi2(is) - ww(ix) * pk2 * j2(as(is)*k) 
     &           * k**3 / (2.d0*pi*pi)
            xi4(is) = xi4(is) + ww(ix) * pk4 * j4(as(is)*k) 
     &           * k**3 / (2.d0*pi*pi)
c            
         enddo
c
      enddo
c
      end
c
c ******************************************************* c
c
      function j0(x)
c
c ******************************************************* c
c
      real*8 j0, x
c     ----------------------------
c
      if(dabs(x).le.1.d-3) then 
         j0 = 1 - x**2/6.d0 + x**4/120.d0
      else
         j0 = sin(x)/x
      endif
c      
      end
c
c ******************************************************* c
c
      function j2(x)
c
c ******************************************************* c
c
      real*8 j2, x
c     ----------------------------
c
      if(dabs(x).le.1.d-3) then 
         j2 = x**2/15.d0 - x**4/210.d0
      else
         j2 = (3.d0-x*x)*dsin(x) - 3.d0*x*dcos(x)
         j2 = j2 / x**3
      endif
c      
      end
c
c ******************************************************* c
c
      function j4(x)
c
c ******************************************************* c
c
      real*8 j4, x
c     ----------------------------
c
      if(dabs(x).le.5.d-1) then 
         j4 = ( 1.d0/945.d0 - x*x/20790.d0 ) * x**4  
      else
         j4 = ( 105.d0 - 45.d0*x*x + x**4 ) * dsin(x) 
     &        - x * ( 105.d0 - 10.d0*x*x ) * dcos(x)
         j4 = j4 / x**5
      endif
c      
      end
c
c ******************************************************* c
c
      subroutine find_pkred(ilabel, kk, pk)
c
c ******************************************************* c
c
c     ilabel =1:  pk0(k)
c     ilabel =2:  pk2(k)
c     ilabel =3:  pk4(k)
c
      implicit none
      integer ilabel, ik_max, ikmax, j, jmin, jmax
      parameter(ikmax=5000)
      real*8  ak(ikmax), pk0(ikmax), pk2(ikmax), pk4(ikmax)
      real*8  kk, s, ds, pk
      common /pkred_data/ ak, pk0, pk2, pk4, ik_max
c     -------------------------------------------
c
      call hunt(ak, ik_max, kk, j)
c
      jmin = j - 2
      jmax = j + 2
      if(jmin.lt.1) jmin = 1
      if(jmax.ge.ik_max) jmax = ik_max
c
      if(ilabel.eq.1) 
     &     call polint(ak(jmin),pk0(jmin),jmax-jmin+1,kk,s,ds)
      if(ilabel.eq.2) 
     &     call polint(ak(jmin),pk2(jmin),jmax-jmin+1,kk,s,ds)
      if(ilabel.eq.3) 
     &     call polint(ak(jmin),pk4(jmin),jmax-jmin+1,kk,s,ds)
c
      pk = s
c      
      end
c
c ******************************************************* c
c
      subroutine pk_red(ik, sigmav2, pk0, pk2, pk4)
c
c ******************************************************* c
c
      implicit none
      integer ik, ikmax, ik_max, imodel_s
      real*8 sigmav2, pk0, pk2, pk4
      real*8 f, ff, alpha, pi
      parameter(ikmax=5000)
      real*8  ak(ikmax), pk_dd(ikmax), pk_dt(ikmax), pk_tt(ikmax)
      real*8  coeff1, coeff3, coeff5, coeff7, coeff9, gammp
      real*8  AA, BB, CC
      real*8  Omega_b, Omega_m, Omega_v, h, Tcmb, n_s, zred
      common  /cosm_param/ Omega_b, Omega_m, h, Tcmb, n_s
      common  /redshift/ zred
      common  /pk_data/ ak, pk_dd, pk_dt, pk_tt, ik_max
      common  /model_s/ imodel_s
c     ----------------------------------------------
c
      Omega_v = 1.d0 - Omega_m
c
      ff = f(zred, Omega_m, Omega_v)
      alpha = ( ak(ik) * ff )**2 * sigmav2
      pi = 4.d0 * datan(1.d0)
c
      if(imodel_s.eq.2) goto  100
c
      if(alpha.ge.0.05d0) then
c
         coeff1 = gammp(0.5d0, alpha) * dsqrt(pi/alpha) 
         coeff3 = gammp(1.5d0, alpha) * dsqrt(pi/alpha**3) * 0.5d0
         coeff5 = gammp(2.5d0, alpha) * dsqrt(pi/alpha**5) * 0.75d0
         coeff7 = gammp(3.5d0, alpha) * dsqrt(pi/alpha**7) * 1.875d0
         coeff9 = gammp(4.5d0, alpha) * dsqrt(pi/alpha**9) * 6.5625d0
c
      elseif(alpha.ge.0.d0 .and. alpha.lt.0.05d0) then
c
         coeff1 = 1.d0 / 0.5d0 - alpha / 1.5d0
         coeff3 = 1.d0 / 1.5d0 - alpha / 2.5d0
         coeff5 = 1.d0 / 2.5d0 - alpha / 3.5d0
         coeff7 = 1.d0 / 3.5d0 - alpha / 4.5d0
         coeff9 = 1.d0 / 4.5d0 - alpha / 5.5d0
c     
      else
c
         write(6,*) ' wrong sign of sigmav2 !!  stop'
         write(6,*) 'alpha, sigmav2=',alpha, sigmav2 
         stop
c
      endif
c
      if(imodel_s.eq.0) then
         pk0 = ( coeff1 * pk_dd(ik) + coeff3 * 2.d0 * ff * pk_dt(ik)
     &        + coeff5 * ff * ff * pk_tt(ik) ) * 0.5d0
         pk2 = ( coeff3 * pk_dd(ik) + coeff5 * 2.d0 * ff * pk_dt(ik)
     &        + coeff7 * ff * ff * pk_tt(ik) ) * 3.75d0 
         pk2 = pk2 - 2.5d0 * pk0
         pk4 = ( coeff5 * pk_dd(ik) + coeff7 * 2.d0 * ff * pk_dt(ik)
     &        + coeff9 * ff * ff * pk_tt(ik) ) * 19.6875d0
         pk4 = pk4 - 4.5d0 * pk2 - 7.875d0 * pk0
      elseif(imodel_s.eq.1) then
         pk0 = ( coeff1 + coeff3 * 2.d0 * ff + coeff5 * ff * ff ) 
     &        * 0.5d0 * pk_dd(ik) 
         pk2 = ( coeff3 + coeff5 * 2.d0 * ff + coeff7 * ff * ff ) 
     &        * 3.75d0 * pk_dd(ik) 
         pk2 = pk2 - 2.5d0 * pk0
         pk4 = ( coeff5 + coeff7 * 2.d0 * ff + coeff9 * ff * ff ) 
     &        * 19.6875d0 * pk_dd(ik) 
         pk4 = pk4 - 4.5d0 * pk2 - 7.875d0 * pk0
      endif
c
 100  if(imodel_s.eq.2) then
c
         if(alpha.le.0.05d0) then 
            AA = 1.d0 - alpha/3.d0 + alpha**2/5.d0
            BB = 1.d0 - 3.d0*alpha/5.d0 + 3.d0*alpha**2/7.d0
            CC = 1.d0 - 5.d0*alpha/7.d0 + 5.d0*alpha**2/9.d0
            pk0 = AA * pk_dd(ik) + 2.d0/3.d0 * ff * BB * pk_dt(ik)
     &           + 1.d0/5.d0 * ff*ff * CC * pk_tt(ik)
c               
            AA = -2.d0*alpha/3.d0 + 4.d0*alpha**2/7.0
            BB = 4.d0/3.d0 - 8.d0*alpha/7.d0 + 20.d0*alpha**2/21.d0
            CC = 4.d0/7.d0 - 10.d0*alpha/21.d0 + 40.d0*alpha**2/99.d0
            pk2 = AA * pk_dd(ik) + ff * BB * pk_dt(ik)
     &           + ff*ff * CC * pk_tt(ik)
            pk4 = pk2           ! this is not entirely correct. 
c
         elseif(alpha.ge.0.05d0) then
            AA = datan(dsqrt(alpha)) / dsqrt(alpha)
            BB = 3.d0/alpha * (1.d0 - AA)
            CC = 5.d0/(3.d0*alpha) * (1.d0 - BB)
c
            pk0 = AA * pk_dd(ik) + 2.d0/3.d0 * ff * BB * pk_dt(ik) 
     &           + 1.d0/5.d0 * ff*ff * CC * pk_tt(ik)
            pk2 = 5.d0/2.d0 * (BB - AA) * pk_dd(ik) 
     &           + ff * ( 4.d0/3.d0 * BB + 3.d0 * ( CC-BB ) ) 
     &           * pk_dt(ik) + ff*ff * ( 3.d0/(2.d0*alpha)*(1.d0-CC) 
     &           - CC/2.d0 ) * pk_tt(ik)
            pk4 = pk2           ! this is not entirely correct. 
         endif
      endif
c
      end
c
c ******************************************************* c
c
      function pk_lin_red(i, pk)
c
c ******************************************************* c
c
c     i = 0:  monopole  
c     i = 2:  quadrupole
c     i = 4:  hexadecapole
c
      implicit none
      integer i
      real*8  pk_lin_red, pk, f, ff
      real*8  Omega_b, Omega_m, Omega_v, h, Tcmb, n_s, zred
      common  /cosm_param/ Omega_b, Omega_m, h, Tcmb, n_s
      common  /redshift/ zred
c     ----------------------------------------------
c
      Omega_v = 1.d0 - Omega_m
c
      ff = f(zred, Omega_m, Omega_v)
c
      if(i.eq.0) pk_lin_red = (1.d0 + 2.d0/3.d0*ff + 1.d0/5.d0*ff*ff )
     &     * pk
      if(i.eq.2) pk_lin_red = ( 4.d0/3.d0 * ff + 4.d0/7.d0 * ff*ff )
     &     * pk
      if(i.eq.4) pk_lin_red = 8.d0/35. * ff*ff 
     &     * pk
c
      end
c
c ******************************************************* c
c
      subroutine calc_pkred(k, a, pkred, dyda, na)
c
c ******************************************************* c
c
      implicit none
      integer ina, na, imodel_s
      real*8  k, kk, sigmav, pkred, dyda(na), a(na)
      real*8  pi, f, ff, alpha, d_sigmav
      parameter(d_sigmav=0.01)
      real*8  zred, Omega_m, Omega_b, h, Tcmb, n_s
      real*8  pk_dd, pk_dt, pk_tt, pkred_up, pkred_down
      real*8  pk_red0, pk_red2, gammp
      real*8  AA, BB, CC
      real*8  coeff1, coeff3, coeff5, coeff7
      common /cosm_param/ Omega_b, Omega_m, h, Tcmb, n_s
      common  /redshift/ zred
      common  /model_s/ imodel_s
c     -------------------------------------------
      pi = 4.d0 * datan(1.d0)
c
      do ina = 1, na
         sigmav = a(ina)
      enddo
c
      ff = f(zred, Omega_m, 1.d0-Omega_m)
c      
      if(k.ge.1.d4) kk = k - 1.d4
      if(k.lt.1.d4) kk = k
c
      call find_pk(1, kk, pk_dd) 
      call find_pk(2, kk, pk_dt) 
      call find_pk(3, kk, pk_tt) 
c
      do ina =1, 3
c
         if(ina.eq.1) alpha = ( kk * ff * sigmav )**2 
         if(ina.eq.2) alpha = ( kk * ff * (sigmav + d_sigmav) )**2 
         if(ina.eq.3) then 
            alpha = ( kk * ff * (sigmav - d_sigmav) )**2 
            if(sigmav.lt.d_sigmav) goto 100
         endif
c
         if(imodel_s.eq.2) goto 10
c
         if(alpha.ge.0.1d0) then
c
            coeff1 = gammp(0.5d0, alpha) * dsqrt(pi/alpha) 
            coeff3 = gammp(1.5d0, alpha) * dsqrt(pi/alpha**3) * 0.5d0
            coeff5 = gammp(2.5d0, alpha) * dsqrt(pi/alpha**5) * 0.75d0
            coeff7 = gammp(3.5d0, alpha) * dsqrt(pi/alpha**7) * 1.875d0
c
         elseif(alpha.ge.0.d0 .and. alpha.lt.0.1d0) then
c
            coeff1 = 1.d0 / 0.5d0 - alpha / 1.5d0
            coeff3 = 1.d0 / 1.5d0 - alpha / 2.5d0
            coeff5 = 1.d0 / 2.5d0 - alpha / 3.5d0
            coeff7 = 1.d0 / 3.5d0 - alpha / 4.5d0
c     
         else
c
            write(6,*) ' wrong sign of sigmav2 !!  stop'
            write(6,*) 'alpha, sigmav2=',alpha, sigmav**2 
            stop
c
         endif
c
         if(imodel_s.eq.0) then
            pk_red0 = ( coeff1 * pk_dd + coeff3 * 2.d0 * ff * pk_dt
     &           + coeff5 * ff * ff * pk_tt ) * 0.5d0
            pk_red2 = ( coeff3 * pk_dd + coeff5 * 2.d0 * ff * pk_dt
     &           + coeff7 * ff * ff * pk_tt ) * 3.75d0 
            pk_red2 = pk_red2 - 2.5d0 * pk_red0
         elseif(imodel_s.eq.1) then
            pk_red0 = ( coeff1 + coeff3 * 2.d0 * ff + coeff5 * ff * ff ) 
     &           * 0.5d0 * pk_dd 
            pk_red2 = ( coeff3 + coeff5 * 2.d0 * ff + coeff7 * ff * ff ) 
     &           * 3.75d0 * pk_dd
            pk_red2 = pk_red2 - 2.5d0 * pk_red0
         endif
c
 10      if(imodel_s.eq.2) then
            if(alpha.le.0.05d0) then
               AA = 1.d0 - alpha/3.d0 + alpha**2/5.d0
               BB = 1.d0 - 3.d0*alpha/5.d0 + 3.d0*alpha**2/7.d0
               CC = 1.d0 - 5.d0*alpha/7.d0 + 5.d0*alpha**2/9.d0
               pk_red0 = AA * pk_dd + 2.d0/3.d0 * ff * BB * pk_dt
     &              + 1.d0/5.d0 * ff*ff * CC * pk_tt
c
               AA = -2.d0*alpha/3.d0 + 4.d0*alpha**2/7.0
               BB = 4.d0/3.d0 - 8.d0*alpha/7.d0 + 20.d0*alpha**2/21.d0
               CC = 4.d0/7.d0 - 10.d0*alpha/21.d0 + 40.d0*alpha**2/99.d0
               pk_red2 = AA * pk_dd + ff * BB * pk_dt
     &              + ff*ff * CC * pk_tt
c     
            elseif(alpha.ge.0.05d0) then
               AA = datan(dsqrt(alpha)) / dsqrt(alpha)
               BB = 3.d0/alpha * (1.d0 - AA)
               CC = 5.d0/(3.d0*alpha) * (1.d0 - BB)
c
               pk_red0 = AA * pk_dd + 2.d0/3.d0 * ff * BB * pk_dt
     &              + 1.d0/5.d0 * ff*ff * CC * pk_tt
               pk_red2 = 5.d0/2.d0 * (BB - AA) * pk_dd
     &              + ff * ( 4.d0/3.d0 * BB + 3.d0 * ( CC-BB ) ) * pk_dt
     &              + ff*ff * ( 3.d0/(2.d0*alpha)*(1.d0-CC) 
     &              - CC/2.d0 ) * pk_tt
            endif
         endif
c
         if(ina.eq.1) then
            if(k.ge.1.d4) pkred = pk_red2
            if(k.lt.1.d4) pkred = pk_red0
         elseif(ina.eq.2) then
            if(k.ge.1.d4) pkred_up = pk_red2
            if(k.lt.1.d4) pkred_up = pk_red0
         elseif(ina.eq.3) then
            if(k.ge.1.d4) pkred_down = pk_red2
            if(k.lt.1.d4) pkred_down = pk_red0
         endif
c
      enddo
c
cccc artificial correction factor so as to match P(k) in linear scale
c
ccc      pkred = 0.976 * pkred
c
 100  do ina = 1, na
         dyda(ina) = (pkred_up - pkred_down ) / d_sigmav / 2.d0
         if(sigmav.lt.d_sigmav) 
     &        dyda(ina) = (pkred_up - pkred ) / d_sigmav 
      enddo
c
      end
c
c ******************************************************* c
c
      subroutine get_initial_sigmav2(sigmav2_guess)
c
c ******************************************************* c
c
      implicit none
      integer ik, ik_max1, ikmax
      parameter(ikmax=5000)
      real*8  sigmav2_guess, pi
      real*8  ak1(ikmax), pk_dd(ikmax), pk_dt(ikmax), pk_tt(ikmax)
      common  /pk_data/ ak1, pk_dd, pk_dt, pk_tt, ik_max1
c     ---------------------------------------------------------
      pi = 4.d0 * datan(1.d0)
c
      sigmav2_guess = 0.d0
c
      do ik = 1, ik_max1
         sigmav2_guess = sigmav2_guess + (pk_tt(ik+1) + pk_tt(ik)) * 
     &        (ak1(ik+1)-ak1(ik)) /2.d0
      enddo
c
      sigmav2_guess = sigmav2_guess / (2.d0*pi*pi) / 3.d0
c
      end
c
c ******************************************************* c
c
      subroutine resorting_nbody_data(ichoice, kmin, kmax, aksim, 
     &     pksim, varsim, ik_max_nbody)
c
c ******************************************************* c
c
      integer ichoice, ikmax, ik_max2, ik_max_nbody
      integer ikk, ik_kmin, ik_kmax 
      parameter(ikmax=5000)
      real*8  ak2(ikmax), pk_sim0(ikmax), pk_sim2(ikmax)
      real*8  var0(ikmax), var2(ikmax) 
      real*8  k, kmin, kmax
      real*8  aksim(ikmax), pksim(ikmax), varsim(ikmax)
      common  /pk_nbody/ ak2, pk_sim0, pk_sim2, var0, var2, ik_max2
c     ------------------------------------------------------------
c
c     find the data point of kmin, kmax
c
      do ik = 1, ik_max2
         if(ak2(ik).le.kmin .and. ak2(ik+1).ge.kmin ) ik_kmin = ik 
         if(ak2(ik).le.kmax .and. ak2(ik+1).ge.kmax ) ik_kmax = ik + 1
      enddo
c
cc      write(6,*) ik_kmin, ik_kmax
cc      write(6,*) ak2(ik_kmin), ak2(ik_kmax)
cc      pause
c
c     re-sorting the data 
c
      ikk = 1
      do ik = 1, ik_max2 
         if(ik.ge.ik_kmin .and. ik.le.ik_kmax) then 
            aksim(ikk) = ak2(ik)
            pksim(ikk) = pk_sim0(ik)
            varsim(ikk) = var0(ik)
ccc            write(6,*) ikk, aksim(ikk), pksim(ikk), varsim(ikk)
            ikk = ikk + 1
         endif
      enddo
c
      ik_max_nbody = ikk -1
c
      if(ichoice.eq.1) then
c
         ikk = ik_max_nbody + 1
         do ik = 1, ik_max2 
            if(ik.ge.ik_kmin .and. ik.le.ik_kmax) then 
               aksim(ikk) = 1.d4 + ak2(ik)
               pksim(ikk) = pk_sim2(ik)
               varsim(ikk) = var2(ik)
ccc               write(6,*) ikk, aksim(ikk), pksim(ikk), varsim(ikk)
               ikk = ikk + 1
            endif
         enddo
c
         ik_max_nbody = ikk - 1
c
      endif
c
      end
c
c ******************************************************* c
c
      subroutine find_pk(ilabel, kk, pk)
c
c ******************************************************* c
c
c     ilabel =1:  Pk_dd(k)
c     ilabel =2:  Pk_dt(k)
c     ilabel =3:  Pk_tt(k)
c
      implicit none
      integer ilabel, ik_max, ikmax, j, jmin, jmax
      parameter(ikmax=5000)
      real*8  ak(ikmax), pk_dd(ikmax), pk_dt(ikmax), pk_tt(ikmax)
      real*8  kk, s, ds, pk
      common  /pk_data/ ak, pk_dd, pk_dt, pk_tt, ik_max
c     -------------------------------------------
c
      call hunt(ak, ik_max, kk, j)
c
      jmin = j - 2
      jmax = j + 2
      if(jmin.lt.1) jmin = 1
      if(jmax.ge.ik_max) jmax = ik_max
c
      if(ilabel.eq.1) 
     &     call polint(ak(jmin),pk_dd(jmin),jmax-jmin+1,kk,s,ds)
      if(ilabel.eq.2) 
     &     call polint(ak(jmin),pk_dt(jmin),jmax-jmin+1,kk,s,ds)
      if(ilabel.eq.3) 
     &     call polint(ak(jmin),pk_tt(jmin),jmax-jmin+1,kk,s,ds)
c
      pk = s
c      
      end
c
c ******************************************************* c
c
      function growth(zred, Omega_m, Omega_v)
c
c ******************************************************* c
c
c     Linear growth factor for flat cosmology
c
      implicit none
ccc      real*8  growth, zred, zred1
ccc      real*8  Omega_m, Omega_v, Omega_k
ccc      real*8  Omega_mz, Omega_vz, gz, g0
      real*8 growth, zred, Omega_m, Omega_v, a, b, c
      real*8 zred1, zi, zz, gz, g0
c     --------------------------------------------
c
ccc      zred1 = 1.d0 + zred
ccc      Omega_k = 1.d0 - Omega_m - Omega_v
c
ccc      Omega_mz = Omega_m * zred1**3 /
ccc     &     ( zred1**2 * Omega_k + zred1**3 * Omega_m + Omega_v )
ccc      Omega_vz = Omega_v /
ccc     &     ( zred1**2 * Omega_k + zred1**3 * Omega_m + Omega_v )
c
ccc      gz = 2.5 * Omega_mz / ( Omega_mz**(4./7.) - Omega_vz + 
ccc     &     (1.d0 + Omega_mz / 2.d0 ) * (1.d0 +Omega_vz / 70.d0 ) )
ccc      g0 = 2.5 * Omega_m / ( Omega_m**(4./7.) - Omega_v + 
ccc     &     (1.d0 + Omega_m / 2.d0 ) * (1.d0 +Omega_v / 70.d0 ) )
c
      a = 1.d0/3.d0
      b = 1.d0
      c = 11.d0/6.d0
c
      zred1 = 1.d0 + zred
      zi = - Omega_v / Omega_m
      zz = zi / zred1**3 
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
      function f(zred, Omega_m, Omega_v)
c
c ******************************************************* c
c
c     d lnD_+(z) / d ln a,   as function of redshift 
c
      implicit none
c
      real*8  f, zred, zred1, Omega_m, Omega_v
      real*8  a, b, c, zi, zz, g1, g2
cc      real*8 eps, growth
c      parameter(eps=1.d-4)
c     ---------------------------------------------------
c     
cc      if( (zred-eps) .ge. 0.d0 ) then 
cc         f = ( dlog(growth(zred+eps, Omega_m, Omega_v)) 
cc     &        - dlog(growth(zred-eps, Omega_m, Omega_v)) ) 
cc     &        / (2.d0 * eps) 
cc         f = - f * (1.d0 + zred)
cc      else
cc         f = ( dlog(growth(zred+eps**2, Omega_m, Omega_v)) 
cc     &        - dlog(growth(zred, Omega_m, Omega_v)) ) / eps**2 
cc         f = - f * (1.d0 + zred)
cc      endif
c
      zred1 = 1.d0 + zred
      zi = - Omega_v / Omega_m
      zz = zi / zred1**3 
c
      a = 4.d0/3.d0
      b = 2.d0
      c = 17.d0/6.d0
      call HYGFX(a,b,c,zz, g1)  
c
      a = 1.d0/3.d0
      b = 1.d0
      c = 11.d0/6.d0
      call HYGFX(a,b,c,zz, g2)  
c
      f = 1.d0 + 6.d0/11.d0 * zz * g1 / g2
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
c
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
c *************************************************************
c
      SUBROUTINE mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,nca,chisq,
     *funcs,alamda)
c
c *************************************************************
c
      implicit none
      INTEGER ma,nca,ndata,ia(ma),MMAX
      REAL*8  alamda,chisq,a(ma),alpha(nca,nca),covar(nca,nca),
     *sig(ndata),x(ndata),y(ndata)
      PARAMETER (MMAX=20)
      external funcs
CU    USES covsrt,gaussj,mrqcof
      INTEGER j,k,l,m,mfit
      REAL*8  ochisq,atry(MMAX),beta(MMAX),da(MMAX)
      SAVE ochisq,atry,beta,da,mfit
      if(alamda.lt.0.)then
        mfit=0
        do 11 j=1,ma
          if (ia(j).ne.0) mfit=mfit+1
11      continue
        alamda=0.001
        call mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nca,chisq,funcs)
        ochisq=chisq
        do 12 j=1,ma
          atry(j)=a(j)
12      continue
      endif
      j=0
      do 14 l=1,ma
        if(ia(l).ne.0) then
          j=j+1
          k=0
          do 13 m=1,ma
            if(ia(m).ne.0) then
              k=k+1
              covar(j,k)=alpha(j,k)
            endif
13        continue
          covar(j,j)=alpha(j,j)*(1.+alamda)
          da(j)=beta(j)
        endif
14    continue
      call gaussj(covar,mfit,nca,da,1,1)
      if(alamda.eq.0.)then
        call covsrt(covar,nca,ma,ia,mfit)
        return
      endif
      j=0
      do 15 l=1,ma
        if(ia(l).ne.0) then
          j=j+1
          atry(l)=a(l)+da(j)
        endif
15    continue
      call mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,nca,chisq,funcs)
      if(chisq.lt.ochisq)then
        alamda=0.1*alamda
        ochisq=chisq
        j=0
        do 17 l=1,ma
          if(ia(l).ne.0) then
            j=j+1
            k=0
            do 16 m=1,ma
              if(ia(m).ne.0) then
                k=k+1
                alpha(j,k)=covar(j,k)
              endif
16          continue
            beta(j)=da(j)
            a(l)=atry(l)
          endif
17      continue
      else
        alamda=10.*alamda
        chisq=ochisq
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software z!0(0.
c
c *************************************************************
c
      SUBROUTINE mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nalp,chisq,
     *funcs)
c
c *************************************************************
c
      implicit none
      INTEGER ma,nalp,ndata,ia(ma),MMAX
      REAL*8 chisq,a(ma),alpha(nalp,nalp),beta(ma),sig(ndata),x(ndata),
     *y(ndata)
      EXTERNAL funcs
      PARAMETER (MMAX=20)
      INTEGER mfit,i,j,k,l,m
      REAL*8 dy,sig2i,wt,ymod,dyda(MMAX)
      mfit=0
      do 11 j=1,ma
        if (ia(j).ne.0) mfit=mfit+1
11    continue
      do 13 j=1,mfit
        do 12 k=1,j
          alpha(j,k)=0.
12      continue
        beta(j)=0.
13    continue
      chisq=0.
      do 16 i=1,ndata
        call funcs(x(i),a,ymod,dyda,ma)
        sig2i=1./(sig(i)*sig(i))
        dy=y(i)-ymod
        j=0
        do 15 l=1,ma
          if(ia(l).ne.0) then
            j=j+1
            wt=dyda(l)*sig2i
            k=0
            do 14 m=1,l
              if(ia(m).ne.0) then
                k=k+1
                alpha(j,k)=alpha(j,k)+wt*dyda(m)
              endif
14          continue
            beta(j)=beta(j)+dy*wt
          endif
15      continue
        chisq=chisq+dy*dy*sig2i
16    continue
      do 18 j=2,mfit
        do 17 k=1,j-1
          alpha(k,j)=alpha(j,k)
17      continue
18    continue
      return
      END
c
c *************************************************************
c
      SUBROUTINE covsrt(covar,npc,ma,ia,mfit)
c
c *************************************************************
c
      implicit none
      INTEGER ma,mfit,npc,ia(ma)
      REAL*8 covar(npc,npc)
      INTEGER i,j,k
      REAL*8 swap
      do 12 i=mfit+1,ma
        do 11 j=1,i
          covar(i,j)=0.
          covar(j,i)=0.
11      continue
12    continue
      k=mfit
      do 15 j=ma,1,-1
        if(ia(j).ne.0)then
          do 13 i=1,ma
            swap=covar(i,k)
            covar(i,k)=covar(i,j)
            covar(i,j)=swap
13        continue
          do 14 i=1,ma
            swap=covar(k,i)
            covar(k,i)=covar(j,i)
            covar(j,i)=swap
14        continue
          k=k-1
        endif
15    continue
      return
      END
c
c
c *************************************************************
c
      SUBROUTINE gaussj(a,n,np,b,m,mp)
c
c
c *************************************************************
c
      implicit none
      INTEGER m,mp,n,np,NMAX
      REAL*8 a(np,np),b(np,mp)
      PARAMETER (NMAX=50)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      REAL*8 big,dum,pivinv
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                pause 'singular matrix in gaussj'
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      END
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
