c-----------------------------------------------------------------------
c © A J S Hamilton 2001
c-----------------------------------------------------------------------
      subroutine gphi(angle,rp,cm,np,rpi,cmi,tol,phi,iord)
      integer np
      real*10 angle,rp(3,np),cm(np),rpi(3),cmi,tol
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
      real*10 angtol,big
c        local variables
      integer i,iseg,j,jm,jml,jmu,jp,jpl,jpu,km,kp,ni,scmi
      real*10 dph,p,ph,phm,php,xi(3),yi(3)
c *
c * Angle along circle about unit direction rpi satisfying
c *    1 - r.rpi = cmi
c * and bounded by
c *    1 - r.rp(j) <= cm(j)  (if cm(j).ge.0)
c *    1 - r.rp(j) > -cm(j)  (if cm(j).lt.0)
c * for j=1,np where rp(j) are unit directions.
c * If the circle lies along a border of the polygon,
c * then the returned angle is zero.
c *
c  Input: rp(3,j),j=1,np
c         cm(j),j=1,np
c         np
c         rpi(3)
c         cmi
c         tol
c Output: angle
c Work arrays: phi and iord should be dimensioned at least 2*np
c
c        set azimuthal angle of non-intersection to big
      data big /1.e6_10/
c        ok if angle tests not too far outside [0,max]
      data angtol /1.e-10_10/
c
c        initialise angle to zero
      angle=0._10
c        check for null circle
      if (cmi.lt.0._10) goto 410
      if (cmi.gt.2._10) goto 410
c        check for zero angle because one circle is null
      if (gzeroar(cm,np).eq.0) goto 410
      scmi=1
c........construct cartesian axes with z-axis along rpi
      call gaxisi(rpi,xi,yi)
c........angles phi about z-axis rp(i) of intersection of i & j circles
c        passing i=0 means circle at edge is considered outside polygon
      call gphij(rp,cm,np,0,rpi,scmi,cmi,xi,yi,big,tol,ni,phi)
c        i circle lies outside (or at edge of) polygon
      if (ni.le.-1) goto 410
c........order angles around circle
c        circle has no intersections
      if (ni.eq.0) then
        angle=TWOPI
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
c         print *,'segment',km,kp,' dph/(2*pi)=',dph/TWOPI
          angle=angle+dph
c        do another segment
        goto 220
      endif
  280 continue
c........check angle is between 0 and 2*pi
      p=angle/TWOPI
c     print *,rpi(1),rpi(2),rpi(3),'angle/(2*pi) =',p
      if (p.lt.0._10) then
        write (*,'(" *** from gphi: angle/(2*pi) = ",g24.16,
     *    " should be >= 0")') p
        goto 420
      elseif (p.gt.1._10) then
c        check if discrepancy is from numerical roundoff
        if (angle.le.TWOPI+angtol) then
          angle=TWOPI
        else
          write (*,'(" *** from gphi: angle/(2*pi) = ",g24.16,
     *      " should be <= 1")') p
          goto 420
        endif
      endif
c........done
      return
c
c        zero angle
  410 continue
      return
c
  420 print *,'*** from gphi: total failure at tol =',tol
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
