c-----------------------------------------------------------------------
c © A J S Hamilton 2001
c-----------------------------------------------------------------------
      subroutine gcmlim(rp,cm,np,rpi,cmimin,cmimax,tol,phi,iord)
      integer np
      real*10 rp(3,np),cm(np),rpi(3),cmimin,cmimax,tol
c        work arrays (could be automatic if compiler supports it)
      integer iord(2*np)
      real*10 phi(2,np)
c
c        parameters
      include 'pi.par'
c        intrinsics
      intrinsic abs
c        externals
      integer gsegij,gzeroar
c        data variables
      real*10 big
c        local variables
      integer i,iseg,jm,jml,jmu,jp,jpl,jpu,ni,scmi
      integer km,kp
      logical inmax,inmin
      real*10 cmi,cmik,cmim,dph,ph,phm,php,phimax,phimin,
     *  si,sik,xi(3),yi(3)
c *
c * Minimum and maximum values of cmi = 1-cos(th)
c * between unit direction rpi
c * and surface of sphere of unit radius bounded by
c *    1 - r.rp(i) < cm(i)  (if cm(i).ge.0)
c *    1 - r.rp(i) > -cm(i)  (if cm(i).lt.0)
c * for i=1,np where rp(i) are unit directions.
c *
c  Input: rp(3,i),i=1,np
c         cm(i),i=1,np
c         np
c         rpi(3)
c         tol
c Output: cmimin, cmimax = min, max values of cmi
c            < 0 means region encloses limiting circle
c            > 0 means region excludes limiting circle
c Work arrays: phi and iord should be dimensioned at least 2*np
c
      data big /1.e6_10/
c
c        check for zero area because one circle is null
      if (gzeroar(cm,np).eq.0) goto 410
      cmimin=2._10
      cmimax=0._10
      inmin=.true.
      inmax=.true.
c--------identify boundary segments around each circle i in turn
      do 280 i=1,np
c        cm(i).ge.2 means include whole sphere, which is no constraint
        if (cm(i).ge.2._10) goto 280
c        scmi * cmi = 1-cos th(i)
        if (cm(i).ge.0._10) then
          scmi=1
        else
          scmi=-1
        endif
        cmi=abs(cm(i))
c        si = sin th(i)
        si=sqrt(cmi*(2._10-cmi))
c        cmik = 1-cos th(ik), th(ik)=angle twixt rpi & rp(i)
        cmik=((rpi(1)-rp(1,i))**2+(rpi(2)-rp(2,i))**2
     *    +(rpi(3)-rp(3,i))**2)/2._10
c        sik = sin th(ik)
        sik=sqrt(cmik*(2._10-cmik))
c        min circle is outside area
        if ((cm(i).ge.0._10.and.cmik.ge.cmi)
     *    .or.(cm(i).lt.0._10.and.cmik.le.cmi)) inmin=.false.
c        max circle is outside area
        if ((cm(i).ge.0._10.and.cmik.le.2._10-cmi)
     *    .or.(cm(i).lt.0._10.and.cmik.ge.2._10-cmi)) inmax=.false.
c........cartesian axes with z-axis along rp(i), x-axis towards rpi
        call gaxisii(rpi,rp(1,i),xi,yi)
c........angles phi about z-axis rp(i) of intersection of i & j circles
        call gphij(rp,cm,np,i,rp(1,i),scmi,cmi,xi,yi,big,tol,ni,phi)
c        i circle lies outside polygon
        if (ni.eq.-1) goto 280
c        area of polygon is zero
        if (ni.eq.-2) goto 410
c........i circle has no intersections
        if (ni.eq.0) then
c        reduce cmimin?
          cmim=cmi+cmik-cmi*cmik-si*sik
          if (cmim.lt.cmimin) cmimin=cmim
c        increase cmimax?
          cmim=cmi+cmik-cmi*cmik+si*sik
          if (cmim.gt.cmimax) cmimax=cmim
c........i circle has intersections
        elseif (ni.gt.0) then
c        find ordering of intersection angles around i circle
          call findbot(phi,2*np,iord,ni)
c        phimin, max are nearest, furthest points from rpi
          phimin=big
          phimax=big
c........vertices around i circle
          jpl=0
c        come here to do another segment
  220     continue
c........is segment edge of polygon?
            iseg=gsegij(rp,cm,np,0,0,i,rp(1,i),scmi,cmi,tol,ni,
     *        phi,iord,jml,jmu,jpl,jpu,1,jm,jp,km,kp,phm,php,ph,dph)
c        error
            if (iseg.eq.-1) goto 420
c        not an edge
            if (iseg.eq.0) goto 220
c        gone full circle
            if (iseg.eq.2) goto 240
            if (php.ge.phm) then
c        segment contains nearest point in i circle, phi=0
              if (phm.le.0._10.and.php.ge.0._10) phimin=0._10
            elseif (php.lt.phm) then
c        segment contains nearest point in i circle, phi=0
              if (phm.le.0._10.or.php.ge.0._10) phimin=0._10
c        segment contains furthest point in i circle, phi=pi
              phimax=0._10
            endif
c        check if segment endpoints tighten limits
            phm=abs(phm)
            php=abs(php)
            if (phm.lt.phimin) phimin=phm
            if (php.lt.phimin) phimin=php
            if (PI-phm.lt.phimax) phimax=PI-phm
            if (PI-php.lt.phimax) phimax=PI-php
c        do another segment
          goto 220
c        reduce cmimin?
  240     if (phimin.ne.big) then
            cmim=cmi+cmik-cmi*cmik-si*sik*cos(phimin)
            if (cmim.lt.cmimin) cmimin=cmim
          endif
c        increase cmimax?
          if (phimax.ne.big) then
            cmim=cmi+cmik-cmi*cmik+si*sik*cos(phimax)
            if (cmim.gt.cmimax) cmimax=cmim
          endif
        endif
  280 continue
c        region encloses limiting circle
      if (inmin) cmimin=-cmimin
      if (inmax) cmimax=-cmimax
      return
c
c        null area
  410 cmimin=2._10
      cmimax=2._10
      return
c
  420 print *,'*** from gmclim: total failure at tol =',tol
      return
c
      end
c
