c ********************************************************** c
c
      program read_pk2
c
c                     Time-stamp: <2011-06-15 19:11:05 ataruya> 
c ********************************************************** c
c
c     Reading the output data of non-linear P(k) from stdPT2.f and/or 
c     pk_closure.f, and estimating non-linear power spectrum or 
c     correlation function in both real and redshift spaces
c
c     Multipole moments of power spectra and correlation function up to 
c     l=10 in redshift space are computed based on the model of 
c     Scoccimarro (2004) (see Eq.(71) of his paper)
c
c     1. Format of output data, pkstd2.dat or pkrenorm2_*.dat : 
c     
c     k, pklin_{no-wiggle}, pklin, P_dd, P_dt, P_tt
c
c     2. Format of output data, pkstd_red.dat or pkrenorm_*_red.dat : 
c
c     k, pklin_{r,no-wiggle}, pklin_r, P_dd, pklin_{0,no-wiggle}, pklin_0, 
c     P0, pklin_{2,no-wiggle}, pklin_2, P2, pklin_{4,no-wiggle}, pklin_4, 
c     P4, P6, P8, P10
c
c     3. Format of output data, xirenorm_*_red.dat (optional): 
c
c     r, xilin_r, Xi_r, xi0_lin, Xi0, xi2_lin, Xi2, xi4_lin, 
c     Xi4, Xi6, Xi8, Xi10
c
      implicit none
c
      integer ik, ikmax, ik_max, itype, idata, iparams, isigmav2
      integer is, is_max, ialpha, imodel
      parameter(ikmax=5000)
      real*8  ak(ikmax), ak2(ikmax), pk_EH(ikmax), pk22_tt(ikmax)
      real*8  pk_lin(ikmax), pk13_dd(ikmax), pk13_dt(ikmax)
      real*8  pk13_tt(ikmax), pk22_dd(ikmax), pk22_dt(ikmax)
      real*8  pk_dd(ikmax), pk_dt(ikmax), pk_tt(ikmax)
      real*8  pk_2loop_dd, pk_2loop_dt, pk_2loop_tt
      real*8  ss, tt, uu, vv, ww, xx, zred, dd2
      real*8  Tcmb, sigma8, n_s, Omega_m, Omega_b, Omega_v, h, w_de
      real*8  growth, f, ff, kmin, kmax
      real*8  sigmav2, pkr_lin(ikmax), pkr(ikmax)
      real*8  pk_red0(ikmax), pk_red2(ikmax), pk_red4(ikmax)
      real*8  pk_red6(ikmax), pk_red8(ikmax), pk_red10(ikmax)
      real*8  pk_lin_red0(ikmax), pk_lin_red2(ikmax), pk_lin_red4(ikmax)
      real*8  pk_EH_red0(ikmax), pk_EH_red2(ikmax), pk_EH_red4(ikmax)
      real*8  dummy, var0(ikmax), var2(ikmax), var4(ikmax), num
      real*8  pk_sim0(ikmax), pk_sim2(ikmax), pk_sim4(ikmax)
      real*8  as(ikmax), xir_lin(ikmax), xir(ikmax)
      real*8  xi0_lin(ikmax), xi2_lin(ikmax), xi4_lin(ikmax)
      real*8  xi0(ikmax), xi2(ikmax), xi4(ikmax)
      real*8  xi6(ikmax), xi8(ikmax), xi10(ikmax)
      character infile*60
      common /cosm_param/ Omega_b, Omega_m, h, Tcmb, n_s
      common /DarkEnergy_EOS/  w_de
      common /pk_tt_data/ ak, pk_lin
      common /pk_data/  ak2, pkr_lin, pkr, pk_lin_red0, pk_red0,
     &     pk_lin_red2, pk_red2, pk_lin_red4, pk_red4, pk_red6,
     &     pk_red8, pk_red10, ik_max
      common /kmin_kmax/  kmin, kmax  
c     ---------------------------------------------------
c     
c     Cosmological parameters for linear growh factor (default values)
c
      h = 0.701d0
      Tcmb = 2.726d0
      n_s = 0.96d0
      sigma8 = 0.817d0
      Omega_m = 0.279d0
      Omega_b = 0.165 * Omega_m
      Omega_v = 1.d0-Omega_m 
      w_de = -1.d0
c
c     ///// Choice of data set ///// c
c
      write(6,*) 'Choose data:                       '
      write(6,*) '                     stdPT2.dat [0]'
      write(6,*) '       renormalized_pk2_RPT.dat [1]'
      write(6,*) '   renormalized_pk2_Closure.dat [2]'
      read(5,*) idata
      if((idata.lt.0) .or. (idata.gt.2) ) stop
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
      write(6,*) 'Note--. spatial curvature is forced to be flat'
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
c ///// Read data file (1) ///// c    
c
 8    if(idata.eq.0) open(9,file='stdPT2.dat',status='unknown')
      if(idata.eq.1) 
     &     open(9,file='renormalized_pk2_RPT.dat',status='unknown')
      if(idata.eq.2) 
     & open(9,file='renormalized_pk2_Closure.dat',status='unknown')
      if(idata.eq.3) then 
         write(6,*) ' input file name ' 
         read(5,*) infile
         open(9,file=infile,status='unknown')
      endif
c
      do ik = 1, ikmax
         read(9,*,END=10) ak(ik), pk_lin(ik), 
     &        pk13_dd(ik), pk13_dt(ik), pk13_tt(ik), 
     &        pk22_dd(ik), pk22_dt(ik), pk22_tt(ik)
         ak2(ik) = ak(ik)
      enddo
c
 10   ik_max = ik - 1
      if(idata.eq.1 .or. idata.eq.2) then
         call calc_sigmav2(ik_max, ss)
         sigmav2 =  ss
      endif
c
      close(9)
c
c ///// Read data file (2) ///// c
c
      if(idata.eq.1 .or. idata.eq.2) then
         write(6,*) 'Do you use alpha-corrected propagator ? y[0],n[1]'
         read(5,*) ialpha
         if(ialpha.eq.0) then
            if(idata.eq.1) 
     &           open(9,file='renormalized_pk2_RPT_alpha.dat',
     &           status='unknown')
            if(idata.eq.2)
     &           open(9,file='renormalized_pk2_Closure_alpha.dat',
     &           status='unknown')
c
            do ik =1, ikmax
               read(9,*,END=15) ak(ik), dummy, 
     &              pk13_dd(ik), pk13_dt(ik), pk13_tt(ik), 
     &              dummy, dummy, dummy
            enddo
 15         ik_max = min(ik - 1, ik_max)
            close(9)
         endif
      endif
c
c ///// Read data file (3) ///// c
c
      if(idata.eq.1 .or. idata.eq.2) then
         write(6,*) ' Do you include 2-loop data ? y[0], n[1] '
         read(5,*) itype
         if(itype.eq.0) then
            if(idata.eq.1) 
     &           open(9,file='renormalized_pk2_RPT_2loop.dat',
     &           status='unknown')
            if(idata.eq.2)
     &           open(9,file='renormalized_pk2_Closure_2loop.dat',
     &           status='unknown')
c
            do ik =1, ikmax
               read(9,*,END=20) ak(ik), pk_lin(ik), 
     &              pk_2loop_dd, pk_2loop_dt, pk_2loop_tt 
               pk22_dd(ik) = pk22_dd(ik) + pk_2loop_dd
               pk22_dt(ik) = pk22_dt(ik) + pk_2loop_dt
               pk22_tt(ik) = pk22_tt(ik) + pk_2loop_tt
            enddo
 20         ik_max = min(ik - 1, ik_max)
            close(9)
         endif
      endif
c
c
      write(6,*) 'ik_max=',ik_max
c      pause
c
c     ///// linear growth rate ///// c
c
      write(6,*) ' type (0, z) or (1, D(z)^2) '
      read(5,*) itype, ss
      if(itype.ne.0 .and.  itype.ne.1) stop
c
      if(itype.eq.0) then 
         zred = ss
         dd2 = growth(zred, Omega_m, Omega_v)**2 
      elseif(itype.eq.1) then 
         dd2 = ss
      endif
c
c     ///// Linear P(k) of no-wiggles : Eisenstein & Hu (1998) ///// c
c
      call no_wiggle_Pk(sigma8, ik_max, ak, pk_EH)
c
c     //////// Save data (1) real space ////// c
c
      if(idata.eq.0) open(10,file='pkstd2.dat',status='unknown')
      if(idata.eq.1) open(10,file='pkrenorm2_RPT.dat',status='unknown')
      if(idata.eq.2) 
     &     open(10,file='pkrenorm2_Closure.dat',status='unknown')
c
      if(idata.eq.0) then
         do ik = 1, ik_max
            pk_dd(ik) = pk_lin(ik) + dd2 * ( pk13_dd(ik) + pk22_dd(ik) )
            pk_dt(ik) = pk_lin(ik) + dd2 * ( pk13_dt(ik) + pk22_dt(ik) )
            pk_tt(ik) = pk_lin(ik) + dd2 * ( pk13_tt(ik) + pk22_tt(ik) )
            pk_dd(ik) = pk_dd(ik) * dd2
            pk_dt(ik) = pk_dt(ik) * dd2
            pk_tt(ik) = pk_tt(ik) * dd2
            pk_lin(ik) = pk_lin(ik) * dd2
            write(10,'(1p6e18.10)') ak(ik), dd2*pk_EH(ik), 
     &           pk_lin(ik), pk_dd(ik), pk_dt(ik), pk_tt(ik)
         enddo
      elseif(idata.eq.1 .or. idata.eq.2) then
         do ik = 1, ik_max
            pk_dd(ik) = pk13_dd(ik) + pk22_dd(ik) 
            pk_dt(ik) = pk13_dt(ik) + pk22_dt(ik) 
            pk_tt(ik) = pk13_tt(ik) + pk22_tt(ik) 
            write(10,'(1p6e18.10)') 
     &           ak(ik), dd2*pk_EH(ik), pk_lin(ik), pk_dd(ik), 
     &           pk_dt(ik), pk_tt(ik) 
         enddo
      endif
c
      close(10)
c
      write(6,*)
      if(idata.eq.0) write(6,*) 'output file --> pkstd2.dat '
      if(idata.eq.1) write(6,*) 'output file --> pkrenorm2_RPT.dat'
      if(idata.eq.2) write(6,*) 'output file --> pkrenorm2_Closure.dat'
      write(6,*)
c
c     //////// Save data (2) redshift space ////// c
c
      if(idata.eq.0) open(11,file='pkstd_red.dat',status='unknown')
      if(idata.eq.1) 
     &     open(11,file='pkrenorm_RPT_red.dat',status='unknown')
      if(idata.eq.2) 
     &     open(11,file='pkrenorm_Closure_red.dat',status='unknown')
      if(idata.eq.3) open(11,file='pknbody_red.dat',status='unknown')
c
c
c     //////// Calculating velocity dispersion ////// c
c
      ff = f(zred, Omega_m, Omega_v)
c
      if(idata.eq.0) then
         call calc_sigmav2(ik_max, ss)
         sigmav2 =  ss
      endif
c
      write(6,*) 'f(z)=', f(zred, Omega_m, Omega_v)
      write(6,*) 'sigmav2=',sigmav2
      write(6,'(A,1p1e12.6,A,1p1e12.6,A)') 
     &        'kmin=',ak(1), ', kmax=',ak(ik_max)
      write(6,*) ' use this value ? y[0], n[1]'
      read(5,*)  isigmav2
      if(isigmav2.gt.1 .or. isigmav2.lt.0) stop
      if (isigmav2.eq.1) then
         write(6,*) 'input  sigmav2'
         read(5,*) sigmav2
      endif
c     
      write(6,*)
      write(6,*) ' Select model of redshift-space distortion:'
      write(6,*) '             Scoccimarro (2004) [0]        '
      write(6,*) '             Kaiser + exp-FoG   [1]        '
      read(5,*)  imodel
c     
      do ik = 1, ik_max
c
         call calc_pk_red(imodel, sigmav2, ak(ik), pk_dd(ik), 
     &        pk_dt(ik), pk_tt(ik), Omega_m, Omega_v, zred, 
     &        ss, tt, uu, vv, ww, xx)
         pkr(ik) = pk_dd(ik)
         pk_red0(ik) = ss
         pk_red2(ik) = tt
         pk_red4(ik) = uu
         pk_red6(ik) = vv
         pk_red8(ik) = ww
         pk_red10(ik) = xx
         pkr_lin(ik) = pk_lin(ik)
         pk_lin_red0(ik) = (1.d0 + 2.d0/3.d0*ff + 1.d0/5.d0*ff*ff )
     &        * pk_lin(ik) 
         pk_lin_red2(ik) = ( 4.d0/3.d0 * ff + 4.d0/7.d0 * ff*ff )
     &        * pk_lin(ik)
         pk_lin_red4(ik) =  8.d0/35. * ff*ff * pk_lin(ik)
         pk_EH_red0(ik) = (1.d0 + 2.d0/3.d0*ff + 1.d0/5.d0*ff*ff )
     &        * pk_EH(ik) * dd2
         pk_EH_red2(ik) = ( 4.d0/3.d0 * ff + 4.d0/7.d0 * ff*ff )
     &        * pk_EH(ik) * dd2
         pk_EH_red4(ik) = 8.d0/35. * ff*ff * pk_EH(ik) * dd2
c
         write(11,'(1p16e18.10)') ak(ik), pk_EH(ik) * dd2, 
     &        pkr_lin(ik), pkr(ik), pk_EH_red0(ik),
     &        pk_lin_red0(ik), pk_red0(ik), pk_EH_red2(ik), 
     &        pk_lin_red2(ik), pk_red2(ik), pk_EH_red4(ik), 
     &        pk_lin_red4(ik), pk_red4(ik), pk_red6(ik), 
     &        pk_red8(ik), pk_red10(ik)
c
      enddo
c
      close(11)
c
      write(6,*)
      if(idata.eq.0) write(6,*) 'output file --> pkstd_red.dat '
      if(idata.eq.1) write(6,*) 'output file --> pkrenorm_RPT_red.dat'
      if(idata.eq.2) write(6,*) 
     &     'output file --> pkrenorm_Closure_red.dat'
c
c     ////////// Output redshift xi(s) data //////////
c
      if(idata.eq.1 .or. idata.eq.2) then
         write(6,*) 
         write(6,*) ' Save also correlation functions ? Yes[0], No[1]'
         read(5,*) iparams
         if(iparams.eq.1) then
            write(6,*) ' Bye ! '
            stop
         endif
c
         kmin = ak(1)
         kmax = ak(ik_max)
c
         call calc_xi(as, xir_lin, xir, xi0_lin, xi0, xi2_lin, xi2, 
     &        xi4_lin, xi4, xi6, xi8, xi10, is_max)
c
         if(idata.eq.1) 
     &        open(10, file='xirenorm_RPT_red.dat',status='unknown')
         if(idata.eq.2) 
     &        open(10, file='xirenorm_Closure_red.dat',status='unknown')
c     
         do is = 1, is_max
            write(10,'(1p12e18.10)') as(is), xir_lin(is), xir(is), 
     &           xi0_lin(is), xi0(is), xi2_lin(is), xi2(is), 
     &           xi4_lin(is), xi4(is), xi6(is), xi8(is), xi10(is)
         enddo
c
         close(10)
c
         write(6,*)
         if(idata.eq.1) write(6,*) 
     &        'output file --> xirenorm_RPT_red.dat'
         if(idata.eq.2) write(6,*) 
     &        'output file --> xirenorm_Closure_red.dat'
c
      endif
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
      real*8 growth, zred, Omega_m, Omega_v, w_de, a, b, c
      real*8 zred1, zi, zz, gz, g0
      common /DarkEnergy_EOS/  w_de
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
c ************************************************ c
c
      subroutine no_wiggle_Pk(sigma8, ik_max, ak, pk_EH)
c
c ************************************************ c
c
      implicit none
      integer ik, ikmax, ik_max
      parameter(ikmax=5000)
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
c
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
      common /cosm_param/ omegab, omega0, h, Tcmb, n_s
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
      subroutine calc_sigmav2(ik_max, sigmav2)
c
c ******************************************************* c
c
      implicit none
      integer ikmax, ik, ik_max
      parameter(ikmax=5000)
      real*8  ak(ikmax), pk(ikmax), sigmav2, pi
      common /pk_tt_data/ ak, pk
c     -------------------------------------------------
      pi = 4.d0 * datan(1.d0)
c
      sigmav2 = 0.d0
c
      do ik=1, ik_max-1
c
         sigmav2 = sigmav2 + (pk(ik+1) + pk(ik)) * 
     &        (ak(ik+1)-ak(ik)) /2.d0
c
      enddo
c
      sigmav2 = sigmav2 /(2*pi*pi) /3.d0
c
      end
c
c ******************************************************* c
c
      subroutine calc_pk_red(imodel, sigmav2, k, pk_dd, pk_dt, pk_tt, 
     &     Omega_m, Omega_v, zred, pk_red0, pk_red2, pk_red4, pk_red6,
     &     pk_red8, pk_red10)
c
c ******************************************************* c
c
c     Model of redshift-space distortion: 
c
c       imodel = 0:  Scoccimarro (2004)
c       imodel = 1:  Kaiser + Exp-FoG
c
      implicit none
      integer imodel
      real*8 sigmav2, k, pk_dd, pk_dt, pk_tt
      real*8 Omega_m, Omega_v, zred, f, ff
      real*8 pk_red0, pk_red2, pk_red4, pk_red6, pk_red8, pk_red10
      real*8 fact, pi, alpha
      common /alpha_param/ alpha
c     -----------------------------------------
c
      pi = 4.d0 * datan(1.d0)
c
      ff = f(zred, Omega_m, Omega_v)
      alpha = ( k * ff )**2 * sigmav2
c
      if(imodel.eq.0) then
         pk_red0 =  fact(0,0)* pk_dd + 2.d0*ff*fact(1,0)* pk_dt
     &        + ff*ff*fact(2,0)*pk_tt
         pk_red2 =  fact(0,2)* pk_dd + 2.d0*ff*fact(1,2)* pk_dt
     &        + ff*ff*fact(2,2)*pk_tt
         pk_red4 =  fact(0,4)* pk_dd + 2.d0*ff*fact(1,4)* pk_dt
     &        + ff*ff*fact(2,4)*pk_tt
         pk_red6 =  fact(0,6)* pk_dd + 2.d0*ff*fact(1,6)* pk_dt
     &        + ff*ff*fact(2,6)*pk_tt
         pk_red8 =  fact(0,8)* pk_dd + 2.d0*ff*fact(1,8)* pk_dt
     &        + ff*ff*fact(2,8)*pk_tt
         pk_red10 =  fact(0,10)* pk_dd + 2.d0*ff*fact(1,10)* pk_dt
     &        + ff*ff*fact(2,10)*pk_tt
      elseif(imodel.eq.1) then
         pk_red0 =  ( fact(0,0) + 2.d0*ff*fact(1,0)
     &        + ff*ff*fact(2,0) )*pk_dd
         pk_red2 =  ( fact(0,2) + 2.d0*ff*fact(1,2)
     &        + ff*ff*fact(2,2) )*pk_dd
         pk_red4 =  ( fact(0,4) + 2.d0*ff*fact(1,4)
     &        + ff*ff*fact(2,4) )*pk_dd
         pk_red6 =  ( fact(0,6) + 2.d0*ff*fact(1,6)
     &        + ff*ff*fact(2,6) )*pk_dd
         pk_red8 =  ( fact(0,8) + 2.d0*ff*fact(1,8)
     &        + ff*ff*fact(2,8) )*pk_dd
         pk_red10 =  ( fact(0,10) + 2.d0*ff*fact(1,10)
     &        + ff*ff*fact(2,10) )*pk_dd
      endif
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
            fact = 12.*gamhalf(n)*gammp(0.5+nn,alpha)/alpha**(nn+0.5)  
     &           -120.*gamhalf(n+1)*gammp(1.5+nn,alpha)/alpha**(nn+1.5)
     &           +140.*gamhalf(n+2)*gammp(2.5+nn,alpha)/alpha**(nn+2.5)
            fact = fact * 9./128.
         elseif(l.eq.6) then
            fact = 
     &          13.*( -5.*gamhalf(n)*gammp(0.5+nn,alpha)/alpha**(nn+0.5)
     &          +105.*gamhalf(n+1)*gammp(1.5+nn,alpha)/alpha**(nn+1.5) )
     &           -4095.*gamhalf(n+2)*gammp(2.5+nn,alpha)/alpha**(nn+2.5)
     &           +3003.*gamhalf(n+3)*gammp(3.5+nn,alpha)/alpha**(nn+3.5)
            fact = fact / 64.
         elseif(l.eq.8) then
            fact = 35.*gamhalf(n)*gammp(0.5+nn,alpha)/alpha**(nn+0.5)
     &           -1260.*gamhalf(n+1)*gammp(1.5+nn,alpha)/alpha**(nn+1.5) 
     &           +6930.*gamhalf(n+2)*gammp(2.5+nn,alpha)/alpha**(nn+2.5)
     &          -12012.*gamhalf(n+3)*gammp(3.5+nn,alpha)/alpha**(nn+3.5)
     &           +6435.*gamhalf(n+4)*gammp(4.5+nn,alpha)/alpha**(nn+4.5)
            fact = fact * 17./512.
         elseif(l.eq.10) then
            fact = -63.*gamhalf(n)*gammp(0.5+nn,alpha)/alpha**(nn+0.5)
     &           +3465.*gamhalf(n+1)*gammp(1.5+nn,alpha)/alpha**(nn+1.5) 
     &          -30030.*gamhalf(n+2)*gammp(2.5+nn,alpha)/alpha**(nn+2.5)
     &          +90090.*gamhalf(n+3)*gammp(3.5+nn,alpha)/alpha**(nn+3.5)
     &         -109395.*gamhalf(n+4)*gammp(4.5+nn,alpha)/alpha**(nn+4.5)
     &          +46189.*gamhalf(n+5)*gammp(5.5+nn,alpha)/alpha**(nn+5.5)
            fact = fact * 21./1024.
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
         elseif(l.eq.6) then
            fact = 52.*dble(n*(2-3*n+n**2))/
     &           dble(105+352*n+344*n**2+128*n**3+16*n**4) -
     &           52.*dble(n*(n**2-1))/
     &           dble(945+1488*n+824*n**2+192*n**3+16*n**4)*alpha +
     &           26.*dble(n*(2+3*n+n*2))/
     &           dble(3465+3776*n+1496*n**2+256*n**3+16*n**4)*alpha**2
         elseif(l.eq.8) then
            fact = 136.*dble(n*(-6+11*n-6*n**2+n**3))/
     &           dble(945+3378*n+3800*n**2+1840*n**3+400*n**4+32*n**5) -
     &           136.*dble(n*(2-n-2*n**2+n**3))/
     &           dble(10395+18258*n+12040*n**2+3760*n**3+560*n**4+
     &           32*n**5)*alpha +
     &           68.*dble(n*(-2-n+2*n**2+n**3))
     &           /dble(45045+56018*n+27000*n**2+6320*n**3+720*n**4+
     &           32*n**5)*alpha**2
         elseif(l.eq.10) then
            fact = 336.*dble(n*(24-50*n+35*n**2-10*n**3+n**4))/
     &           dble(10395+39048*n+48556*n**2+27840*n**3+8080*n**4+
     &           1152*n**5+64*n**6) -
     &           336.*dble(n*(-6+5*n+5*n**2-5*n**3+n**4))/
     &           dble(135135+258144*n+193036*n**2+72960*n**3+14800*n**4+
     &           1536*n**5+64*n**6)*alpha +
     &           168.*dble(n*(4-5*n+n**4))
     &           /dble(675675+930360*n+517036*n**2+148800*n**3+
     &           23440*n**4+1920*n**5+64*n**6)*alpha**2
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
c     Gamma(n+1/2) up to n=8
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
      if(n.eq.9) gamhalf = 67303.564453125d0
      if(n.eq.10) gamhalf = 639383.8623046875d0
      if(n.eq.11) gamhalf = 6.71353055419921875d6
      if(n.eq.12) gamhalf = 7.720560137329101563d7
c
      gamhalf = gamhalf * dsqrt(pi)
c
      end
c
c ******************************************************* c
c
      subroutine calc_xi(as, xir_lin, xir, xi0_lin, xi0, xi2_lin, xi2, 
     &     xi4_lin, xi4, xi6, xi8, xi10, is_max)
c
c ******************************************************* c
c
      implicit none
      integer  ikmax, is, ismax, is_max, ix, ixmax
      parameter(ikmax=5000, ismax=501, ixmax=1000)
      real*8  as(ikmax), xir_lin(ikmax), xir(ikmax)
      real*8  xi0_lin(ikmax), xi2_lin(ikmax), xi4_lin(ikmax)
      real*8  xi0(ikmax), xi2(ikmax), xi4(ikmax)
      real*8  xi6(ikmax), xi8(ikmax), xi10(ikmax)
      real*8  smin, smax, pi, j0, j2, j4, j6, j8, j10, kmin, kmax
      real*8  pkr, pk0, pk2, pk4, pk6, pk8, pk10 
      real*8  pkr_lin, pk0_lin, pk2_lin, pk4_lin
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
         xir_lin(is) = 0.d0
         xi0_lin(is) = 0.d0
         xi2_lin(is) = 0.d0
         xi4_lin(is) = 0.d0
         xir(is) = 0.d0
         xi0(is) = 0.d0
         xi2(is) = 0.d0
         xi4(is) = 0.d0
         xi6(is) = 0.d0
         xi8(is) = 0.d0
         xi10(is) = 0.d0
c
         do ix=1, ixmax
c
            k = dexp(kk(ix)) 
            call find_pk(1, k, pkr_lin)
            call find_pk(2, k, pk0_lin)
            call find_pk(3, k, pk2_lin)
            call find_pk(4, k, pk4_lin)
            call find_pk(5, k, pkr)
            call find_pk(6, k, pk0)
            call find_pk(7, k, pk2)
            call find_pk(8, k, pk4)
            call find_pk(9, k, pk6)
            call find_pk(10, k, pk8)
            call find_pk(11, k, pk10)
c
            xir_lin(is) = xir_lin(is) + ww(ix) * pkr_lin * j0(as(is)*k) 
     &           * k**3 / (2.d0*pi*pi)
            xi0_lin(is) = xi0_lin(is) + ww(ix) * pk0_lin * j0(as(is)*k) 
     &           * k**3 / (2.d0*pi*pi)
            xi2_lin(is) = xi2_lin(is) - ww(ix) * pk2_lin * j2(as(is)*k) 
     &           * k**3 / (2.d0*pi*pi)
            xi4_lin(is) = xi4_lin(is) + ww(ix) * pk4_lin * j4(as(is)*k) 
     &           * k**3 / (2.d0*pi*pi)
            xir(is) = xir(is) + ww(ix) * pkr * j0(as(is)*k) 
     &           * k**3 / (2.d0*pi*pi)
            xi0(is) = xi0(is) + ww(ix) * pk0 * j0(as(is)*k) 
     &           * k**3 / (2.d0*pi*pi)
            xi2(is) = xi2(is) - ww(ix) * pk2 * j2(as(is)*k) 
     &           * k**3 / (2.d0*pi*pi)
            xi4(is) = xi4(is) + ww(ix) * pk4 * j4(as(is)*k) 
     &           * k**3 / (2.d0*pi*pi)
            xi6(is) = xi6(is) + ww(ix) * pk6 * j6(as(is)*k) 
     &           * k**3 / (2.d0*pi*pi)
            xi8(is) = xi8(is) + ww(ix) * pk8 * j8(as(is)*k) 
     &           * k**3 / (2.d0*pi*pi)
            xi10(is) = xi10(is) + ww(ix) * pk10 * j10(as(is)*k) 
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
      function j6(x)
c
c ******************************************************* c
c
      real*8 j6, x
c     ----------------------------
c
      if(dabs(x).le.5.d-1) then 
         j6 =  ( 1.d0/135135.d0 - x*x/4054050.d0 ) * x**6  
      else
         j6 = ( 10395.d0 - 4725.d0*x*x + 210.d0*x**4 -x**6 ) * dsin(x) 
     &        - x * 21.d0*( 495.d0 - 60.d0*x*x + x**4 ) * dcos(x)
         j6 = j6 / x**7
      endif
c      
      end
c
c ******************************************************* c
c
      function j8(x)
c
c ******************************************************* c
c
      real*8 j8, x
c     ----------------------------
c
      if(dabs(x).le.5.d-1) then 
         j8 =  ( 1.d0/34459425.d0 - x*x/1309458150.d0 ) * x**8  
      else

         j8 = (2027025.d0 - 945945.d0*x**2 + 51975.d0*x**4 - 
     &        630.d0*x**6 + x**8)*dsin(x) + 9.d0*x* (-225225.d0 + 
     &        30030.d0*x**2 - 770.d0*x**4 + 4.d0*x**6)*dcos(x)
         j8 = j8 / x**9
      endif
c      
      end
c
c
c ******************************************************* c
c
      function j10(x)
c
c ******************************************************* c
c
      real*8 j10, x
c     ----------------------------
c
      if(dabs(x).le.5.d-1) then 
         j10 =  ( 1.d0/13749310575.d0 - x*x/632468286450.d0 ) * x**10  
      else

         j10 =  (-55.d0*(11904165.d0 - 1670760.d0*x**2 + 51597.d0*x**4 - 
     &        468.d0*x**6 + x**8)*dcos(x))/x**10 + 
     &        ((654729075.d0 - 310134825.d0*x**2 + 18918900.d0*x**4 - 
     &        315315.d0*x**6 + 1485.d0*x**8 - x**10)*dsin(x))/x**11
      endif
c      
      end
c
c ******************************************************* c
c
      subroutine find_pk(ipk, kk, pk)
c
c ******************************************************* c
c
      implicit none
      integer ik_max, ikmax, ipk
      integer j, jmin, jmax
      parameter(ikmax=5000)
      real*8 ak(ikmax), pkr_lin(ikmax), pkr(ikmax)
      real*8 pk_lin_red0(ikmax), pk_lin_red2(ikmax), pk_lin_red4(ikmax)
      real*8 pk_red0(ikmax), pk_red2(ikmax), pk_red4(ikmax)
      real*8 pk_red6(ikmax), pk_red8(ikmax), pk_red10(ikmax)
      real*8 kk, s, ds, pk
      common /pk_data/  ak, pkr_lin, pkr, pk_lin_red0, pk_red0,
     &     pk_lin_red2, pk_red2, pk_lin_red4, pk_red4, pk_red6,
     &     pk_red8, pk_red10, ik_max
c     -------------------------------------------
c
      call hunt(ak, ik_max, kk, j)
c
      jmin = j - 2
      jmax = j + 2
      if(jmin.lt.1) jmin = 1
      if(jmax.ge.ik_max) jmax = ik_max
c
      if(ipk.eq.1) 
     &     call polint(ak(jmin),pkr_lin(jmin),jmax-jmin+1,kk,s,ds)
      if(ipk.eq.2) 
     &     call polint(ak(jmin),pk_lin_red0(jmin),jmax-jmin+1,kk,s,ds)
      if(ipk.eq.3) 
     &     call polint(ak(jmin),pk_lin_red2(jmin),jmax-jmin+1,kk,s,ds)
      if(ipk.eq.4) 
     &     call polint(ak(jmin),pk_lin_red4(jmin),jmax-jmin+1,kk,s,ds)
      if(ipk.eq.5) 
     &     call polint(ak(jmin),pkr(jmin),jmax-jmin+1,kk,s,ds)
      if(ipk.eq.6) 
     &     call polint(ak(jmin),pk_red0(jmin),jmax-jmin+1,kk,s,ds)
      if(ipk.eq.7) 
     &     call polint(ak(jmin),pk_red2(jmin),jmax-jmin+1,kk,s,ds)
      if(ipk.eq.8) 
     &     call polint(ak(jmin),pk_red4(jmin),jmax-jmin+1,kk,s,ds)
      if(ipk.eq.9) 
     &     call polint(ak(jmin),pk_red6(jmin),jmax-jmin+1,kk,s,ds)
      if(ipk.eq.10) 
     &     call polint(ak(jmin),pk_red8(jmin),jmax-jmin+1,kk,s,ds)
      if(ipk.eq.11) 
     &     call polint(ak(jmin),pk_red10(jmin),jmax-jmin+1,kk,s,ds)
c
      pk = s
c      
      end
c
c ******************************************************* c
c
      function f(zred, Omega_m, Omega_v)
c
c ******************************************************* c
c
c     d lnD_+(z) / d ln a,   as function of redshift (exact)
c
      implicit none
c
      real*8  f, zred, zred1, Omega_m, Omega_v, w_de
      real*8  a, b, c, zi, zz, g1, g2
      common /DarkEnergy_EOS/  w_de
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
c
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
