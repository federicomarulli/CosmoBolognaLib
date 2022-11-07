c-----------------------------------------------------------------------
c © A J S Hamilton 2001
c-----------------------------------------------------------------------
      subroutine gvphi(angle,v,rp,cm,np,rpi,cmi,vi,tol,phi,iord)
      integer np
      real*10 angle,v(3),rp(3,np),cm(np),rpi(3),cmi,vi(3),tol
c        work arrays (could be automatic if compiler supports it)
      integer iord(2*np)
      real*10 phi(2,np)
c
c        parameters
      include 'pi.par'
      real*10 TWOPI
      parameter (TWOPI=2._10*PI)
c        externals
      integer gsegij,gzeroar
c        data variables
      real*10 big
c        local (automatic) variables
      integer i,iphin,iseg,j,jm,jml,jmu,jp,jpl,jpu,km,kp,ni,scmi
      real*10 dph,dphin,dphinmn,ph,phin,phm,php,
     *  si,xi(3),xv,yi(3),yv,zv
c *
c * This routine is mostly lifted from gphi and gvert.
c *
c * Point at centre of segment of circle
c *    1 - r.rpi = cmi
c * containing, or otherwise closest to, point vi,
c * lying inside
c *    1 - r.rp(j) <= cm(j)  (if cm(j).ge.0)
c *    1 - r.rp(j) > -cm(j)  (if cm(j).lt.0)
c * for j=1,np where rp(j) are unit directions.
c *
c  Input: rp(3,j),j=1,np
c         cm(j),j=1,np
c         np
c         rpi(3)
c         cmi
c         vi(3) = unit vector desired to lie inside,
c                 or closest to, segment.
c         tol
c Output: angle = angular length of segment of circle
c               = 0. if boundary lies entirely outside circle.
c         v(3) = unit vector at centre of segment of circle.
c Work arrays: phi and iord should be dimensioned at least 2*np
c
      data big /1.e6_10/
c
c        initialise length of segment to zero
      angle=0._10
c        initialise point at centre of segment to zero
      v(1)=0._10
      v(2)=0._10
      v(3)=0._10
c        check for null circle
      if (cmi.lt.0._10) goto 410
      if (cmi.gt.2._10) goto 410
c        check for zero angle because one circle is null
      if (gzeroar(cm,np).eq.0) goto 410
c        initialise dphinmn to impossibly large value
      dphinmn=big
c........si = sin thi
      scmi=1
      si=sqrt(cmi*(2._10-cmi))
c........construct cartesian axes with z-axis along rpi
      call gaxisi(rpi,xi,yi)
c........azimuthal angle closest to vector vi
      xv=xi(1)*vi(1)+xi(2)*vi(2)+xi(3)*vi(3)
      yv=yi(1)*vi(1)+yi(2)*vi(2)+yi(3)*vi(3)
      phin=atan2(yv,xv)
c........angles phi about z-axis rp(i) of intersection of i & j circles
c        passing i=0 means circle at edge is considered outside polygon
      call gphij(rp,cm,np,0,rpi,scmi,cmi,xi,yi,big,tol,ni,phi)
c        i circle lies outside (or at edge of) polygon
      if (ni.le.-1) goto 410
c........order angles around circle
c        circle has no intersections
      if (ni.eq.0) then
        angle=TWOPI
        ph=phin
        xv=si*cos(ph)
        yv=si*sin(ph)
        zv=1._10-cmi
        v(1)=zv*rpi(1)+xv*xi(1)+yv*yi(1)
        v(2)=zv*rpi(2)+xv*xi(2)+yv*yi(2)
        v(3)=zv*rpi(3)+xv*xi(3)+yv*yi(3)
c        circle has intersections
      elseif (ni.gt.0) then
c        find ordering of intersection angles around circle
        call findbot(phi,2*np,iord,ni)
c........vertices around i circle
        jpl=0
c        come here to do another segment
  220   continue
c........is segment edge of polygon?
          iseg=gsegij(rp,cm,np,0,0,i,rpi,scmi,cmi,tol,ni,
     *      phi,iord,jml,jmu,jpl,jpu,1,jm,jp,km,kp,phm,php,ph,dph)
c        error
          if (iseg.eq.-1) goto 420
c        not an edge
          if (iseg.eq.0) goto 220
c        gone full circle
          if (iseg.eq.2) goto 280
c........segment satisfies conditions
c         print *,'segment',km,kp,
c    *      ' (',jm,' in',jml,jmu,',',jp,' in',jpl,jpu,' of',ni,')',
c    *      ' ph=',real(ph),' dph=',real(dph)
c        phase phin to central point ph
          iphin=nint((phin-ph)/TWOPI)
          phin=phin-iphin*TWOPI
          if (phm.le.phin.and.phin.le.php) then
            dphin=0._10
          elseif (phm.gt.phin) then
            dphin=phm-phin
          elseif (phin.gt.php) then
            dphin=phin-php
          endif
c        segment contains or is closest to phin
          if (dphin.lt.dphinmn) then
            dphinmn=dphin
            angle=dph
c        coords of centre of edge in frame where axes are xi, yi, rp
            xv=si*cos(ph)
            yv=si*sin(ph)
            zv=1._10-cmi
            v(1)=zv*rpi(1)+xv*xi(1)+yv*yi(1)
            v(2)=zv*rpi(2)+xv*xi(2)+yv*yi(2)
            v(3)=zv*rpi(3)+xv*xi(3)+yv*yi(3)
c        segment contains phin, so cannot be beaten
            if (dphin.eq.0._10) goto 280
          endif
c        do another segment
        goto 220
      endif
  280 continue
      return
c
c        zero angle
  410 continue
      return
c
  420 print *,'*** from gvphi: total failure at tol =',tol
      write (*,'(a3,a20,4a24)')
     *  ' ','x','y','z',
c    *  'r',
     *  '1-c'
      do j=1,np
        write (*,'(i3,5g24.16)')
     *    j,(rp(i,j),i=1,3),
c    *    sqrt(rp(1,j)**2+rp(2,j)**2+rp(3,j)**2),
     *    cm(j)
      enddo
      write (*,'(i3,5g24.16)')
     *  0,(rpi(i),i=1,3),
c    *  sqrt(rpi(1)**2+rpi(2)**2+rpi(3)**2),
     *  cmi
      return
c
      end
c
