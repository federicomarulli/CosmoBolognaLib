c-----------------------------------------------------------------------
c © A J S Hamilton 2001
c-----------------------------------------------------------------------
      subroutine gsphera(area,bound,vert,w,lmax1,im,nw,ibv,
     *  azmin,azmax,elmin,elmax,v,dw)
      integer lmax1,im,nw,ibv
      real*10 area,bound(2),vert(2),w(im,nw),
     *  azmin,azmax,elmin,elmax,dw(nw)
c        work array (could be automatic if compiler supports it)
      real*10 v(lmax1)
c
c        parameters
      include 'pi.par'
      real*10 TWOPI,PIBYTWO
      parameter (TWOPI=2._10*PI,PIBYTWO=PI/2._10)
c        data variables
      real*10 elmino,elmaxo
c        saved variables
      real*10 cl,cu,dth,sl,su
      save cl,cu,dth,sl,su
c        local (automatic) variables
      integer i,l,m,lm,lmax,mmax
      real*10 azmx,cmph,d,dph,ph,smph,thmin,thmax
c *
c * Accelerated computation of spherical transform
c * of rectangle bounded by lines of constant latitude & longitude.
c *
c  Input: lmax1 = lmax+1 where lmax is maximum desired l of transform.
c         im = 1 means compute only real part of harmonics;
c              2 means compute both real and imaginary parts.
c              Note harmonics are real if region possesses reflection
c              symmetry through plane defined by x- and z-axes.
c         nw = [(lmax+1)*(lmax+2)]/2
c         ibv = 0 to 3 controls evaluation of bound & vert,
c               and determines the sign of the returned harmonics,
c               same as in gspher;
c               see comments in gspher for more details.
c         azmin, azmax = minimum, maximum azimuth of rectangle in radians;
c                if azmin > azmax, it is assumed that the rectangle
c                runs from azmin to azmax + 2*pi .
c                To cover an entire strip of latitude, use
c                azmin, azmax = 0, 2*pi .
c         elmin, elmax = minimum, maximum elevation of rectangle in radians;
c                South pole is at elevation -pi/2, North at +pi/2.
c Output: area = area of rectangle in steradians.
c         bound = length of boundary of rectangle in radians if ibv=0,
c                 or as explained in gspher ibv>0.
c         vert = sum over vertices of 1-psi/tan(psi) if ibv=0,
c                where psi is exterior angle (=pi-interior angle)
c                at vertex, or as explained in gspher if ibv>0.
c         dw = integral_thmin^thmax Y_lm(th,0) sin th d th,
c              which should be saved between calls.
c Input/Output: w(i,lm) = spherical transform, dimensioned w(im,nw)
c            w(i,lm), i=1,im, lm=l*(l+1)/2+m+1, l=0,lmax, m=0,l;
c            w(1,lm) is real part, w(2,lm) is imaginary part (if im=2).
c            Note w(l,-m)=(-)**m*[Complex conjugate of w(l,m)], just as
c                 Y(l,-m)=(-)**m*[Complex conjugate of Y(l,m)].
c Work arrays: v should be dimensioned at least lmax1.
c
      data elmino,elmaxo /2*0._10/
c
c        zero stuff
      area=0._10
      bound(1)=0._10
      bound(2)=0._10
      vert(1)=0._10
      vert(2)=0._10
      do lm=1,nw
        do i=1,im
          w(i,lm)=0._10
        enddo
      enddo
c        check input parameters OK
      if (lmax1.le.0) goto 200
      if (elmin.ge.PIBYTWO.or.elmax.le.-PIBYTWO) goto 200
      if (elmin.ge.elmax) goto 200
      azmx=azmax
c        assume azmax.lt.azmin means need to add 2*pi to azmax
      if (azmx.lt.azmin) azmx=azmx+TWOPI
c--------compute integrals of harmonics if elmin and elmax changed
      if (elmino.ne.elmin.or.elmaxo.ne.elmax) then
        if (elmax.ge.PIBYTWO) then
          thmin=0._10
          cu=1._10
          su=0._10
        else
          thmin=PIBYTWO-elmax
          cu=cos(thmin)
          su=sin(thmin)
        endif
        if (elmin.le.-PIBYTWO) then
          thmax=PI
          cl=-1._10
          sl=0._10
        else
          thmax=PIBYTWO-elmin
          cl=cos(thmax)
          sl=sin(thmax)
        endif
        dth=thmax-thmin
c        integrals of harmonics: this takes most time
        call iylm(thmin,thmax,dw,lmax1,nw,v)
        elmaxo=elmax
        elmino=elmin
      endif
c--------fast computation of harmonics
      dph=azmx-azmin
      area=(cu-cl)*dph
      if (ibv.eq.0.or.ibv.eq.2.or.ibv.eq.3) then
        bound(1)=(sl+su)*dph+2._10*dth
        bound(2)=(1._10/sl-2._10*sl+1._10/su-2._10*su)*dph-2._10*dth
        vert(1)=4._10
        vert(2)=4._10*(cl/sl-cu/su)
        if (ibv.eq.2) then
          bound(1)=-bound(1)
          bound(2)=-bound(2)
          vert(1)=-vert(1)
          vert(2)=-vert(2)
        endif
      endif
      ph=(azmx+azmin)/2._10
      if (ibv.ge.2) dph=-dph
      lmax=lmax1-1
      mmax=lmax
      if (dph-nint(dph/TWOPI)*TWOPI.eq.0._10) mmax=0
      do m=0,mmax
        if (m.eq.0) then
          d=dph
        elseif (m.gt.0) then
          d=sin(m*dph/2._10)*2._10/dble(m)
        endif
        cmph=cos(m*ph)
        smph=sin(m*ph)
        lm=(m*(m+1))/2+1
        do l=m,lmax
          lm=lm+l
          w(1,lm)=w(1,lm)+cmph*d*dw(lm)
          if (im.eq.2) w(2,lm)=w(2,lm)-smph*d*dw(lm)
        enddo
      enddo
  200 continue
      return
      end
c
