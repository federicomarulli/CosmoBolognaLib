c-----------------------------------------------------------------------
c © A J S Hamilton 2001
c-----------------------------------------------------------------------
      real*10 function wrho(az,el,w,lmax,mmax,im,nw,lsmooth,esmooth)
      integer lmax,mmax,im,nw
      real*10 el,az,w(im,nw),lsmooth,esmooth
c
c        parameters
      real*10 HALF
      parameter (HALF=1._10/2._10)
      include 'pi.par'
c        local variables
c     integer i
      integer l,lm,lmax1,l1,m,mmax1,m1
      real*10 al,al1,am,cel,cm,dwrho,lsmoot1,phi,
     *  sel,sm,smooth,t,tm,tp,y,z,zm,zn,zp
c *
c * Given window harmonics w_lm, returns value of window function
c *    sum w_lm Y_lm exp{-[l(l+1)/lsmooth(lsmooth+1)]**(esmooth/2)]}
c * at azimuthal angle az radians, elevation el radians.
c *
c  Input: az = azimuthal angle (= longitude) in radians.
c         el = elevation (= latitude = pi/2 - polar angle) in radians.
c         w(i,lm) = spherical transform, dimensioned w(im,nw)
c            w(i,lm), i=1,im, lm=l*(l+1)/2+m+1, l=0,lmax, m=0,l;
c            w(1,lm) is real part, w(2,lm) is imaginary part (if im>=2).
c         lmax = maximum harmonic number of harmonics to include.
c         mmax = include azimuthal harmonics -mmax to + mmax
c              = lmax normally, but other values are acceptable.
c         im = row dimension of harmonics w.
c         nw = [(lmax+1)*(lmax+2)]/2
c         lsmooth = smoothing harmonic number
c                 = 0 (or < 0) to skip smoothing.
c         esmooth = exponent of smoothing, ignored if lsmooth <= 0
c                 = 2 for Gaussian smoothing.
c Output: wrho = sum w_lm Y_lm
c                    * exp{-[l(l+1)/lsmooth(lsmooth+1)]**(esmooth/2)]}
c
      lmax1=lmax+1
      mmax1=mmax+1
      lsmoot1=lsmooth+1._10
      cel=cos(el)
      sel=sin(el)
      phi=az
      wrho=0._10
      do 180 m1=1,mmax1
        m=m1-1
        am=m
        if (m.eq.0) then
c        zn=z(0,0)
          zn=1._10/sqrt(4._10*PI)
          cm=1._10
          sm=0._10
        elseif (m.gt.0) then
c        zn=z(m,m)
          zn=-sqrt((am-HALF)/am)*cel*zn
          cm=2._10*cos(am*phi)
          if (im.eq.2) sm=2._10*sin(am*phi)
        endif
        z=0._10
        zp=zn
        lm=(m*m1)/2+1
        do 160 l1=m1,lmax1
          al1=l1
          l=l1-1
          al=l
c        zm=z(l-1,m); z=z(l,m); zp=z(l+1,m)
          zm=z
          z=zp
          t=al1+al
          tm=sqrt((al+am)*(al-am))
          tp=sqrt((al1+am)*(al1-am))
          zp=(t*sel*z-tm*zm)/tp
          y=z*sqrt(t)
c        dwrho=w(l,m)*Y(l,m)+w(l,-m)*Y(l,-m)
          lm=lm+l
          dwrho=cm*w(1,lm)
          if (im.eq.2) dwrho=dwrho-sm*w(2,lm)
          dwrho=dwrho*y
          if (lsmooth.gt.0._10) then
            smooth=exp(-(al/lsmooth*al1/lsmoot1)**(esmooth/2._10))
            dwrho=dwrho*smooth
          endif
          wrho=wrho+dwrho
c         write (*,'(2i4,i6,5g12.4)') l,m,lm,sm,cm,y,(w(i,lm),i=1,im)
  160   continue
  180 continue
      return
      end
c
