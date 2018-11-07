
!! Implements the multi-point propagator expansion for the nonlinear       !!
!! matter power spectrum, following Crocce, Scoccimarro and Bernardeau     !! 
!! 2012 (arXiv ). Version 1.0                                              !!
     
      program MPTbreeze
      
      implicit  none

      real*8    tol,pi,ak1,ak2,adk
      integer   iwr,Nbins
      parameter (pi=3.1415927d0)

      real*8    Om0,sig80,sig8,ns,w,var8
      real*8    az,Dp,HF,HF0,sigv,eps,akfun,kc
      real*8    ak,G
      common    / argk    / ak
      common    / klimits / eps,kc
      common    / sigv    / sigv

      integer   i,nlin,nlinmx
      parameter (nlinmx=2000)
      real*8    klin(nlinmx),plin(nlinmx),plind(nlinmx),k,tf,akmin,akmax
      common    / pklin    / klin,plin,plind,nlin
      common    / pklimits / akmin,akmax
      character fileTF*200

      real*8    variance,powlin,af,fk,p1loop,p2loop
      external  variance,powlin,gabqx,af,p1loop,p2loop

      character filePk*200
      real*8    P0,P1L,P2L

      integer   iverbose,nargs,iarg,nc,ic
      real*8    omegab,hv,omch2,ombh2,omc,omb,Hubblev,atf(10)
      character path*60,cambF*200,argc*200
      character args(3)*200,argstr*20


!}{ Some default parameters

      nc       = 2
      iverbose = 1
      path     = "./"
      fileTF   = ""
      filePk   = ""

      sig80    = 1.d4
      az       = 1.d4

!}{ Default cosmological parameters

      Hubblev = 70.
      hv      = 0.7
      sig80   = 0.817
      Om0     = 0.279
      az      = 0.0
      ns      = 0.96
      w       = -1.0

!}{ Reading of the command line options

      nargs   = iargc()

      do iarg=1,nargs
         call getarg(iarg,argc)
         if (argc.eq."-noverbose") then
            iverbose=0
         endif
         if (argc.eq."-verbose") then
            if (iarg.eq.nargs) goto 106
            call getarg(iarg+1,argstr)
            read(argstr,*,err=106) iverbose
         endif
         if (argc.eq."-path") then
            if (iarg.eq.nargs) goto 106
            call getarg(iarg+1,path)
         endif
         if (argc.eq."-fileTF") then
            if (iarg.eq.nargs) goto 106
            call getarg(iarg+1,fileTF)
         endif
         if (argc.eq."-filePk") then
            if (iarg.eq.nargs) goto 106
            call getarg(iarg+1,filePk)
         endif
         if (argc.eq."-camb") then
            if (iarg.eq.nargs) goto 106
            call getarg(iarg+1,cambF)
            cambF=path(1:len_trim(path))//cambF(1:len_trim(cambF))
            ombh2=0.
            omch2=0.
            open(9, file=cambF, status='old',err=104)
            do while(.true.)
              read(9,*,END=14) args(1),args(2),args(3)
              if (args(1).eq.'hubble') read(args(3),*) Hubblev
              if (args(1).eq.'ombh2') read(args(3),*) ombh2
              if (args(1).eq.'omch2') read(args(3),*) omch2
              if (args(1).eq.'omega_baryon') read(args(3),*) omb
              if (args(1).eq.'omega_cdm') read(args(3),*) omc
              if (args(1).eq.'w') read(args(3),*) w
              if (args(1)(1:21).eq.'scalar_spectral_index') 
     &             read(args(3),*) ns
              if (args(1)(1:17).eq.'transfer_filename') then
                 fileTF=args(3)
              endif
              nc=7
            enddo
 14         continue
            close(9)
            hv=Hubblev/100.
            omegab=max(ombh2/hv/hv,omb)
            Om0=omegab+max(omch2/hv/hv,omc)
         endif
         if (argc.eq."-sigma8") then
            if (iarg.eq.nargs) goto 106
            call getarg(iarg+1,argstr)
            read(argstr,*,err=106) sig80
         endif
         if (argc.eq."-omegam") then
            if (iarg.eq.nargs) goto 106
            call getarg(iarg+1,argstr)
            read(argstr,*,err=106) Om0
         endif
         if (argc.eq."-h") then
            if (iarg.eq.nargs) goto 106
            call getarg(iarg+1,argstr)
            read(argstr,*,err=106) hv
         endif
         if (argc.eq."-ns") then
            if (iarg.eq.nargs) goto 106
            call getarg(iarg+1,argstr)
            read(argstr,*,err=106) ns
         endif
         if (argc.eq."-redshift") then
            if (iarg.eq.nargs) goto 106
            call getarg(iarg+1,argstr)
            read(argstr,*,err=106) az
         endif
         if (argc.eq."-w") then
            if (iarg.eq.nargs) goto 106
            call getarg(iarg+1,argstr)
            read(argstr,*,err=106) w
         endif
      enddo

      if(filePk.eq."") filePk="MPTbreeze_"//fileTF(1:len_trim(fileTF))
      fileTF = path(1:len_trim(path))//fileTF(1:len_trim(fileTF))
      filePk = path(1:len_trim(path))//filePk(1:len_trim(filePk))

!}{ Inline instructions when verbose <> 0      

      if (iverbose.ge.1) then
      write(6,*)'-----------------------------------------------------
     &--------------------------'
      write(6,*)':  MPTbreeze  = density power spectrum from Multipoin
     &t Propagators            :'
      write(6,*)':               in Perturbation Theory               
     &                         :'
      write(6,*)'-----------------------------------------------------
     &--------------------------'
      write(6,*)': Description of the online options:                 
     &                         :' 
      write(6,*)':  "-noverbose", to suppress verbose                 
     &                         :'
      write(6,*)':                                                    
     &                         :'
      write(6,*)':  "-path" aaa, access path to input and ouput files 
     & (default =./)           :'
      write(6,*)':  "-camb" aaa, to use CAMB ".ini" file for cosmologi
     &cal parameters           :'
      write(6,*)':           and transfer function file.              
     &                         :'
      write(6,*)':                                                    
     &                         :'
      write(6,*)':  "-fileTF" aaa, name of the input Transfer Function
     &file.                    :'
      write(6,*)':  "-sigma8" xxx, to specify the value of sigma8 (def
     &ault.=0.817)             :'
      write(6,*)':  "-omegam" xxx, to specify the value of Omega_m (de
     &fault=.279)              :'
      write(6,*)':  "-ns" xxx, to specify the value of ns (default=.96
     &)                        :'
      write(6,*)':  "-w" xxx, to specify the value of w (default=-1)  
     &                         :'
      write(6,*)':  "-redshift" z, value of the redshift of the comput
     &ed spectra (default=0)   :'
      write(6,*)':  "-filePk" aaa, name of the ouput power spectra.   
     &                         :'
      write(6,*)':____________________________________________________
     &_________________________:'
      write(6,*)
      write(6,*) '    Execution log :'
      write(6,*)
      endif
      if (iverbose.ge.1) then
         write(6,'(A,A)')     '> Initial Transfer file = ',fileTF(1:len_
     &   trim(fileTF))
         write(6,'(A,f10.4)') '> sigma8 = ',sig80
         write(6,'(A,f10.4)') '> omegam = ',Om0
         write(6,'(A,f10.4)') '> spectral index = ',ns
         write(6,'(A,f10.4)') '> eq. of state of DE = ',w
         write(6,'(A,f10.4)') '> output redshift = ',az
      endif         

!}{ Set up linear power spectrum

      open(unit=11,file=fileTF,status='old',err=105) ! Read Transfer Function
      nlin = 0
      do while (.true.)
         read(11,*,end=13) k,(atf(ic),ic=1,nc-1)
         tf         = atf(nc-1)
         nlin       = nlin+1
         klin(nlin) = k
         plin(nlin) = k**ns*tf**2.d0         
      enddo
 13   continue
      close(11)
      If (iverbose.ge.1) then
         write(6,*)
         write(6,"(A,I5,A20)") '> read ',nlin,' lines from input TF'
      endif

      call spline(klin,plin,nlin,3d30,3d30,plind)
      akmin = klin(1)
      akmax = klin(nlin)
      If (iverbose.ge.1) then
       write(6,"(10x,A,F10.7,A,F8.2)") 'kmin=',akmin,' and kmax=',akmax
      endif


      call GFw(w,0.d0,Om0,HF0)
      call GFw(w,az,Om0,HF)
      Dp = HF/HF0/(1+az)
      sig8 = sig80*Dp
      var8 = variance(8.d0)
      do i=1,nlin
      plin(i) = sig8**2*plin(i)/var8
      plind(i)= sig8**2*plind(i)/var8
      enddo

!}{ Compute BAO damping and nonlinear density power spectrum 

      tol  = 1.d-4                                                    
      iwr  = 0          
      call gabqx(powlin,eps,akmax,sigv,tol,iwr)
      sigv = dsqrt(4.d0*pi*sigv/3.d0)

      adk = 0.005d0
!!    ak1 = 0.01d0
      ak2 = 1.d0/sigv      
      eps = akmin

      If (iverbose.ge.1) then
         write(6,'(10x,A,f7.3,A,f10.4)') 'sigv=',sigv,' and knl=',ak2
      endif

      open(unit=8,file=filePk,status='unknown')      

      i = 1
      do while(ak.le.ak2)         !! Edit here for different maximum wavemode and binning !!
         ak   = adk*real(i)
         kc   = max(20.d0*ak,pi)
         call gabqx(af,eps,kc,fk,tol,iwr)
         G    = dexp(fk)
         write(8,100) ak,powlin(ak)*G**2.d0,p1loop(ak)*G**2.d0,
     &   p2loop(ak)*G**2.d0, fk
         i = i+1
      enddo
      close(8)     
 
      If (iverbose.ge.1) then
      write(6,*)
      write(6,'(A)') '> done '
      write(6,'(A,A)') '> Output PS file = ',filePk(1:len_trim(filePk))
      endif
      
 100  format(2x,8e16.6)
      stop

      goto 107
 104  write(6,*) ' ! Error, file not found: ',cambF(1:len_trim(cambF))
      stop
 105  write(6,*) ' ! Error, file not found: ',fileTF(1:len_trim(fileTF))
      stop
 106  write(6,*) ' ! Error in command line for ',argc(1:len_trim(argc))
      stop
 107  continue

      end program MPTbreeze


![1]

      real*8    function powlin(QQ)   !linear power

      implicit  none
      integer   nlinmx,nlin
      parameter (nlinmx=2000)
      real*8    klin(nlinmx),plin(nlinmx),plind(nlinmx)
      real*8    akmin,akmax,QQ
      common    / pklin / klin,plin,plind,nlin
      common    / pklimits / akmin,akmax
      external  splint
      
      if (QQ.ge.akmin .and. QQ.le.akmax) then 
      call splint(klin,plin,plind,nlin,QQ,powlin)
      else
      powlin=0.d0
      endif

      return
      end

![2]

      real*8    function variance(qq)   !linear variance

      implicit  none
      integer   limit
      parameter (limit=1000)
      integer   neval,ier,iord(limit),last
      real*8    alist(limit),blist(limit),elist(limit),rlist(limit)
      real*8    ans,epsabs,epsrel,abserr
      real*8    akmin,akmax,varint,ar,qq
      common    / varscale / ar
      common    / pklimits / akmin,akmax
      external  varint

      ar = qq
      epsabs = 0.d0
      epsrel = 1.d-4                                                            
      call dqage(varint,akmin,akmax,epsabs,epsrel,30,limit,ans,abserr,
     $neval,ier,alist,blist,rlist,elist,iord,last)
      variance = ans
      return
      end

      real*8    function varint(ak)   !linear variance integrand

      implicit  none
      real*8    ak,ax,wth,pid,ar,powlin
      common    / varscale / ar
      external  powlin
      pid = 3.141592654d0
      ax = ak*ar
      wth = 3.d0/ax**3 *(dsin(ax)-ax*dcos(ax))
      varint = powlin(ak) * 4.d0*pid *ak**2 *wth**2
      return
      end

![3]

      real*8    function af(q)  !f integrand

      implicit  none
      real*8    k,q,pi,powlin
      parameter (pi=3.1415926d0)
      external  powlin
      common    / argk /k
       
      af=(4.96031746031746d-4*(4.d0*(6.d0*k**7*q-7.9d1*k**5*q**
     &3+5.d1*k**3*q**5-2.1d1*k*q**7)+3.d0*(k**2-q**2)**3*(2.d0*
     &k**2+7.d0*q**2)*log((k-q)**2/(k+q)**2)))/(k**3*q**5)
     
      af=af*powlin(q)*4.d0*pi*q**2
      return 
      end

![4]

      real*8    function P1loop(k)   

      implicit  none
      real*8    pi,ans,k,kk,P1loopq,eps,kc
      external  P1loopq,G,gabqx

      common    / klimits / eps,kc
      common    / varsP /   kk

      integer   limit
      parameter (limit=1000)
      integer   neval,ier,iord(limit),last
      real*8    alist(limit),blist(limit),elist(limit),rlist(limit)
      real*8    epsabs,epsrel,abserr

      integer   iwr
      real*8    tol

      pi     = 3.141592654d0
      kk     = k

      epsabs = 0.d0
      epsrel = 0.001d0         !! Relative Error !!

      call dqage(P1loopq,eps,kc,epsabs,epsrel,30,limit,
     $ ans,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      P1loop = 4.d0*pi*ans

      return
      end
      
      real*8    function P1loopq(q)   !  22 q-integrand

      implicit  none
      integer   iwr
      real*8    powlin,qq,tol,P1loopx,ans,q
      common    / loopP / qq
      external  powlin,gabqy,P1loopx

      qq   = q
      tol  = 1.d-2        !! Relative Error !!                  
      iwr  = 0                                                             
      call gabqy(P1loopx,-1.d0,1.d0,ans,tol,iwr)
      P1loopq = q**2*powlin(q)*ans      
      return
      end

      real*8    function P1loopx(x)   !  22 x-integrand

      implicit  none
      real*8    powlin,F2,k1q,a1q,q,x,k
      external  powlin,F2
      common    / varsP / k
      common    / loopP / q

      k1q = dsqrt(k**2+q**2-2.d0*k*q*x)
      a1q = (k*x-q)/k1q
      P1loopx=powlin(k1q)*F2(k1q,q,a1q)**2 
      
      return
      end

![5]

      real*8 function P2loop(k)   ! 2L power spectrum

      include   'header-vegas.f'
      real*8    pi,k,ak,eps,kc
      external  vegas,integrand
      common    / wavevec / ak
      common    / klimits / eps,kc
      
      pi = 3.141592654d0
      ak = k

      call vegas(ndim, ncomp, integrand,
     &  epsrel, epsabs, verbose, mineval, maxeval,
     &  nstart, nincrease,
     &  neval, fail, integral, error, prob)
      P2loop = integral(1)
      
      return
      end

      subroutine integrand(ndim,xx,ncomp,ff)

      implicit  none
      integer   ndim, ncomp
      real*8    xx(ndim),ff(ncomp),f
      real*8    ak,q,p,x,y,phi
      real*8    F3,F3s,powlin,pi,cosxy,kqp,a1,a2,a3
      external  F3,F3s,powlin

      real*8    qmin,qmax,pmin,pmax,xmin,xmax,ymin,ymax,phimin,phimax
      real*8    eps,kc
      common    / klimits / eps,kc
      common    /wavevec/ ak
      
      pi     = 3.141592654d0

      qmin   = eps
      qmax   = kc
      pmin   = eps
      pmax   = kc
      xmin   = -1.d0
      xmax   = 1.d0
      ymin   = -1.d0
      ymax   = 1.d0
      phimin = 0.d0
      phimax = 2.d0*pi

      q   = qmax*xx(1)+qmin*(1.d0-xx(1))
      p   = pmax*xx(2)+pmin*(1.d0-xx(2))
      x   = xmax*xx(3)+xmin*(1.d0-xx(3))
      y   = ymax*xx(4)+ymin*(1.d0-xx(4))
      phi = phimax*xx(5)+phimin*(1.d0-xx(5))

      cosxy = (x*y+dsqrt(1.d0-x**2)*dsqrt(1.d0-y**2)*cos(phi))
      kqp   = dsqrt(ak**2+q**2+p**2-2.*(ak*p*y+ak*q*x-p*q*cosxy))
      
      a1 = (ak*q*x-q**2-p*q*cosxy)/(kqp*q)
      a2 = cosxy
      a3 = (ak*p*y-p**2-p*q*cosxy)/(kqp*p)

      f=F3s(kqp,q,p,a1,a2,a3)**2*powlin(kqp)*powlin(q)*powlin(p)

      f=f*q**2*p**2  

      ff(1)=f*12.d0*pi*(qmax-qmin)*(xmax-xmin)
     &*(pmax-pmin)*(ymax-ymin)*(phimax-phimin)

      return
      end

      real*8  function F3s(q1,q2,q3,x12,x23,x31)   ! symmetrized F3 kernel

      implicit  none
      real*8    q1,q2,q3,x12,x23,x31,F3
      external  F3      
      
      F3s = F3(q1,q2,q3,x12,x23,x31)+F3(q2,q1,q3,x12,x31,x23)+
     &      F3(q3,q2,q1,x23,x12,x31)+F3(q1,q3,q2,x31,x23,x12)+
     &      F3(q2,q3,q1,x23,x31,x12)+F3(q3,q1,q2,x31,x12,x23)
      F3s = F3s/6.d0
      
      return
      end

      real*8 function F3(q1,q2,q3,x12,x23,x31)   ! F3 (needs k=q123 as common)

      implicit none
      real*8   q1,q2,q3,x12,x23,x31,al1,al2,be1,be2,ak,F2,G2
      common   / wavevec / ak
      external F2,G2
      
      al1 = 1.d0+(q2*x12+q3*x31)/q1
      al2 = 1.d0+(q3*q1*x31+q2*q3*x23)/(q1**2+q2**2+2.d0*q1*q2*x12)
      be1 = 0.5d0*ak**2 *(q1*q2*x12+q3*q1*x31)/q1**2
      be1 = be1/(q2**2+q3**2+2.d0*q2*q3*x23)
      be2 = 0.5d0*ak**2 *(q2*q3*x23+q3*q1*x31)/q3**2
      be2 = be2/(q1**2+q2**2+2.d0*q1*q2*x12)
      
      F3  = 7.d0*al1*F2(q2,q3,x23) + 2.d0*be1* G2(q2,q3,x23)+
     &      (7.d0*al2 + 2.d0*be2)*G2(q1,q2,x12) 
      
      F3 =  F3/18.d0

      if(F3.ne.F3) F3 = 0.d0

      return
      end

![6]

      real*8    function F2(q1,q2,x12)   ! F2 kernel

      implicit  none
      real*8    q1,q2,x12
      real*8    epsilon,Om

      F2 = (5.d0/7.d0 + x12/2.d0 *(q1/q2+q2/q1) + 2.d0/7.d0 *x12**2)

      return
      end

      real*8 function G2(q1,q2,x12)      ! G2 kernel

      implicit  none
      real*8    q1,q2,x12

      G2 = 3.d0/7.d0 + x12/2.d0*(q1/q2+q2/q1) + 4.d0/7.d0*x12**2.d0

      return
      end


      include 'subs/dqage.f'
      include 'subs/d1mach.f'
      include 'subs/spline_dp.f'
      include 'subs/splint_dp.f'
      include 'subs/dqawfe.f'      
      include 'subs/gabqx.f'
      include 'subs/gabqy.f'
      include 'subs/gabqz.f'
      include 'subs/GFw.f'
