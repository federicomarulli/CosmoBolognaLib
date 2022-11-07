c-----------------------------------------------------------------------
c © A J S Hamilton 2001
c-----------------------------------------------------------------------
      integer function gzeroar(cm,np)
      integer np
      real*10 cm(np)
c
c        local (automatic) variables
      integer i
c *
c * Check for zero area because one circle is null.
c *
c  Input: cm
c         np
c Return value: 0 if area is zero because one circle is null,
c               1 otherwise
c
      do i=1,np
        if (cm(i).eq.0._10) goto 200
        if (cm(i).le.-2._10) goto 200
      enddo
      gzeroar=1
      return
c
  200 gzeroar=0
      return
c
      end
c
c-----------------------------------------------------------------------
      subroutine gaxisi(rp,xi,yi)
      real*10 rp(3),xi(3),yi(3)
c
c        local (automatic) variables
      real*10 sx
c *
c * Cartesian axes with z-axis along rp.
c *
c  Input: rp
c Output: xi, yi forming right-handed orthonormal system with rp
c
      sx=rp(1)**2+rp(3)**2
      if (sx.gt..5_10) then
        sx=sqrt(sx)
c        xi in direction y x rp (= x direction if rp is along z)
        xi(1)=rp(3)/sx
        xi(2)=0._10
        xi(3)=-rp(1)/sx
      else
        sx=sqrt(rp(1)**2+rp(2)**2)
c        xi in direction rp x z
        xi(1)=rp(2)/sx
        xi(2)=-rp(1)/sx
        xi(3)=0._10
      endif
c        yi in direction rp x xi (= y direction if rp is along z)
      yi(1)=xi(3)*rp(2)-xi(2)*rp(3)
      yi(2)=xi(1)*rp(3)-xi(3)*rp(1)
      yi(3)=xi(2)*rp(1)-xi(1)*rp(2)
      return
c
      end
c
c-----------------------------------------------------------------------
      subroutine gaxisii(rpi,rp,xi,yi)
      real*10 rpi(3),rp(3),xi(3),yi(3)
c
c        local (automatic) variables
      real*10 ri,sik
c *
c * Cartesian axes with z-axis along rp, x-axis towards rpi.
c *
c  Input: rpi
c Output: xi, yi forming right-handed orthonormal system with rpi
c
c        set yi in direction rp x rpi
      yi(1)=rpi(3)*rp(2)-rpi(2)*rp(3)
      yi(2)=rpi(1)*rp(3)-rpi(3)*rp(1)
      yi(3)=rpi(2)*rp(1)-rpi(1)*rp(2)
c        project yi orthogonal to rp, to avoid problems near yi=0
      ri=yi(1)*rp(1)+yi(2)*rp(2)+yi(3)*rp(3)
      yi(1)=yi(1)-ri*rp(1)
      yi(2)=yi(2)-ri*rp(2)
      yi(3)=yi(3)-ri*rp(3)
      sik=yi(1)**2+yi(2)**2+yi(3)**2
      if (sik.gt.0._10) then
c        sik = sin th(ik)
        sik=sqrt(sik)
        yi(1)=yi(1)/sik
        yi(2)=yi(2)/sik
        yi(3)=yi(3)/sik
c        rpi is same/opposite direction to rp: set yi along z x rp
      else
        ri=sqrt(rp(1)**2+rp(2)**2)
        if (ri.gt.0._10) then
          yi(1)=-rp(2)/ri
          yi(2)=rp(1)/ri
          yi(3)=0._10
c        if rp is also along z-axis, set yi along y-axis
        else
          yi(1)=0._10
          yi(2)=1._10
          yi(3)=0._10
        endif
      endif
c        xi in direction yi x rp
      xi(1)=yi(2)*rp(3)-yi(3)*rp(2)
      xi(2)=yi(3)*rp(1)-yi(1)*rp(3)
      xi(3)=yi(1)*rp(2)-yi(2)*rp(1)
      return
c
      end
c
c-----------------------------------------------------------------------
      subroutine gphij(rp,cm,np,i,rpi,scmi,cmi,xi,yi,big,tol,ni,phi)
      integer np,i,scmi,ni
      real*10 rp(3,np),cm(np),rpi(3),cmi,xi(3),yi(3),big,tol,phi(2,np)
c
c        intrinsics
      intrinsic abs
c        local (automatic) variables
      integer j
      real*10 bi,bj,cmij,cmj,d,dc,xj,yj
c *
c * angles phi about z-axis rp(i) of intersection of i & j circles
c * phi = big means no intersection
c *
c  Input: rp, cm
c         np = number of caps
c         i
c         rpi = rp(i)
c         scmi = sign(cm(i))
c         cmi = abs(cm(i))
c         xi, yi form orthonormal axes with rp(i)
c         big
c         tol = great circle angle closer than which
c               i & j circles are considered coincident
c Output: ni = number of intersections of i circle with j circles
c            = -1 if circle i is entirely outside another circle
c            = -2 if area of polygon is zero
c         phi(1,j), phi(2,j) = azimuthal angle about rpi
c              of intersection of j circle with i circle;
c              zero azimuthal angle is in direction xi.
c
c--------initialise phi to big, meaning no intersection
      do j=1,np
        phi(1,j)=big
        phi(2,j)=big
      enddo
c        set count of number of intersections with i circle to zero
      ni=0
c--------find intersection of i circle with each j circle in turn
      do 150 j=1,np
c        skip self
        if (j.eq.i) goto 150
c        cm(j).ge.2 means include whole sphere, so no intersection
        if (cm(j).ge.2._10) goto 150
c        cmij = 2 sin^2[th(ij)/2] = 1-cos th(ij)
        cmij=((rpi(1)-rp(1,j))**2+(rpi(2)-rp(2,j))**2
     *    +(rpi(3)-rp(3,j))**2)/2._10
c        cmj = 1-cos th(j)
c        bj = cj-ci*cij
c        d = 1-ci^2-cj^2-cij^2+2*ci*cj*cij
c        dph = atan(sqrt(d)/bj) is angle from rp(j) to intersection
        cmj=abs(cm(j))
        bj=(cmi-cmj)+cmij*(1._10-cmi)
        d=-(cmi-cmj)**2+cmij*(2._10*((cmi+cmj)-cmi*cmj)-cmij)
c        if i and j circles are angle e apart at closest approach, then
c        d approx 2 sin th(i) sin th(j) sin th(ij) * e for small e
        dc=2._10*sqrt(cmi*cmj*cmij*(2._10-cmi)*(2._10-cmj)*(2._10-cmij))
c........positive d means i and j circles intersect
c        if i and j circles are <= tol apart, treat d as zero
        if (d.gt.tol*dc) then
          d=sqrt(d)
c        ph = atan(yj/xj) is angle from xi to rp(j)
          xj=xi(1)*rp(1,j)+xi(2)*rp(2,j)+xi(3)*rp(3,j)
          yj=yi(1)*rp(1,j)+yi(2)*rp(2,j)+yi(3)*rp(3,j)
c        order intersection angles so segment from first to second angle
c        is inside j circle
c Notice order of evaluation of RHS of phi(1,j) and phi(2,j) is same;
c this ensures that gcc evaluates identically for identical arguments.
          if (cm(j).ge.0._10) then
c        phi(1,j)=ph-dph , phi(2,j)=ph+dph
            phi(1,j)=atan2(yj*bj-xj*d,xj*bj+yj*d)
            phi(2,j)=atan2(yj*bj+xj*d,xj*bj-yj*d)
          elseif (cm(j).lt.0._10) then
c        phi(1,j)=ph+dph , phi(2,j)=ph-dph
            phi(2,j)=atan2(yj*bj-xj*d,xj*bj+yj*d)
            phi(1,j)=atan2(yj*bj+xj*d,xj*bj-yj*d)
          endif
c        increment count of number of intersections
          ni=ni+2
c........zero d means i and j circles just touch;
c        negative d means i and j circles don't intersect
        else
c        bi = ci-cj*cij
          bi=(cmj-cmi)+cmij*(1._10-cmj)
c. . . . bi=0 means i and j circles coincide, implying also bj=0 and d=0
c        but test both bi and bj to guard against numerics
          if (bi.eq.0._10.or.bj.eq.0._10) then
c        null intersection of areas:
c        rp(i) and rp(j) point in same direction
            if (cmij.lt.1._10) then
c        cm(i) and cm(j) have opposite sign
              if ((scmi.ge.0.and.cm(j).lt.0._10)
     *          .or.(scmi.lt.0.and.cm(j).ge.0._10)) goto 220
c        rp(i) and rp(j) point in opposite directions
            elseif (cmij.gt.1._10) then
c        cm(i) and cm(j) have same sign
              if ((scmi.ge.0.and.cm(j).ge.0._10)
     *          .or.(scmi.lt.0.and.cm(j).lt.0._10)) goto 220
            endif
c        only do later of the two degenerate circles
            if (i.lt.j) goto 210
c. . . . i circle does not coincide with j circle
          else
c. . . . i circle is outside j circle
            if ((cm(j).ge.0._10.and.bj.gt.0._10)
     *        .or.(cm(j).lt.0._10.and.bj.lt.0._10)) then
c        j circle also outside i circle means null intersection area
              if ((scmi.ge.0.and.bi.gt.0._10)
     *          .or.(scmi.lt.0.and.bi.lt.0._10)) goto 220
c        skip i circle since it's entirely outside j circle
              goto 210
            endif
c. . . . i circle is inside j circle, and just touches it
            if (d.ge.-tol*dc) then
c. . . . j circle is also inside i circle, and just touches it
              if ((scmi.ge.0.and.bi.lt.0._10)
     *          .or.(scmi.lt.0.and.bi.gt.0._10)) then
c        ph = atan(yj/xj) is angle from xi to rp(j)
                xj=xi(1)*rp(1,j)+xi(2)*rp(2,j)+xi(3)*rp(3,j)
                yj=yi(1)*rp(1,j)+yi(2)*rp(2,j)+yi(3)*rp(3,j)
c        phi(1,j)=phi(2,j)=ph
                phi(1,j)=atan2(yj*bj,xj*bj)
                phi(2,j)=phi(1,j)
c        increment count of number of intersections
                ni=ni+2
              endif
            endif
          endif
        endif
  150 continue
c--------normal return
      return
c
c--------circle lies entirely outside another circle
  210 ni=-1
      return
c
c--------area is zero
  220 ni=-2
      return
c
      end
c
c-----------------------------------------------------------------------
      subroutine ggpij(np,gp,i,big,phi)
      integer np,gp(np),i
      real*10 big,phi(2,np)
c
c        local (automatic) variables
      integer j
c *
c * Called after gphij().
c * Identify friends:
c * if circle j intersects circle i, then i and j are friends.
c *
c  Input: np = number of caps.
c         i = circle
c         big = value of phi if i and j do not intersect.
c         phi(1,j), phi(2,j) = azimuthal angle about rp(i)
c              of intersection of j circle with i circle.
c Input/Output:
c         gp(j),j=1,np = which group of friends circle j belongs to:
c               i and j are friends
c               if i and j circles intersect (anywhere).
c
      do j=1,np
c      j circle intersects i circle
        if (phi(1,j).ne.big) then
c      make j group same as i group
          if (gp(i).lt.gp(j)) then
            gp(j)=gp(i)
c      make i group same as j group
          elseif (gp(j).lt.gp(i)) then
            gp(i)=gp(j)
          endif
        endif
      enddo
      return
      end
c
c-----------------------------------------------------------------------
      subroutine ggp(np,gp)
      integer np,gp(np)
c
c        local (automatic) variables
      integer j
c *
c * Called after all friends have been identified.
c * Consolidate groups by letting friends of friends be friends.
c *
c  Input: np = number of caps.
c Input/Output:
c         gp(j),j=1,np = which group of friends circle j belongs to:
c               i and j are friends
c               if i and j circles intersect (anywhere),
c               and friends of friends are friends.
c
      do j=1,np
  200   if (gp(j).ne.gp(gp(j))) then
          gp(j)=gp(gp(j))
          goto 200
        endif
      enddo
      return
      end
c
c-----------------------------------------------------------------------
      integer function gsegij(rp,cm,np,npb,npc,i,rpi,scmi,cmi,tol,ni,
     *  phi,iord,jml,jmu,jpl,jpu,nphbv,jm,jp,km,kp,phm,php,ph,dph)
      integer np,npb,npc,i,scmi,ni,iord(2*np),
     *  jml,jmu,jpl,jpu,nphbv,jm(nphbv),jp(nphbv),km(nphbv),kp(nphbv)
      real*10 rp(3,np),cm(np),rpi(3),cmi,tol,phi(2,np),phm,php,ph,dph
c
c        parameters
      include 'pi.par'
      real*10 TWOPI
      parameter (TWOPI=2._10*PI)
c        intrinsics
      intrinsic abs
c        local (automatic) variables
      integer iphbv,j,jj,k,kk
      logical ismax
      real*10 bik,cmik,cmk,d,psi,psim(2),dphc,si
c        local variables to be saved
      integer jl,ju
      save jl,ju
c *
c * Determine whether next segment of i circle
c * is an edge of the polygon.
c * The thing that complicates this subroutine is the need
c * to deal with near multiple intersections, where several
c * j circles intersect the i circle at almost the same point.
c * Intersections closer than angle tol (great circle separation)
c * are considered coincident.
c *
c * If the segment is an edge, then km and kp are the circles
c * crossing at the lower and upper limits of the edge,
c * with multiple intersections correctly taken into account.
c * The values of phm and php are those corresponding to the
c * intersection with the km and kp circles.
c *
c * Normally nphbv = 1, but nphbv = 2 when called from gphbv.
c * If nphbv = 2, then the intersections of the i circle
c * are collected into two sets, those from the same polygon
c * as the i circle, and those from the abutting polygon.
c *
c  Input: rp
c         cm
c         npb
c         npc
c         i
c         rpi = rp(i)
c         scmi = sign(cm(i))
c         cmi = abs(cm(i))
c         np
c         tol = great circle angle closer than which
c               points on i circle are considered coincident
c         ni = number of intersections of j circles with i circle
c         phi = angles of j circles about i circle
c         iord = order of j circles about i circle
c         jpl = 0 on first call
c             = per output of previous call subsequently
c Output: jml, jmu = points jml to jmu are at lower point
c         jpl, jpu = points jpl to jpu are at upper point
c                    jml <= jmu < jpl <= jpu
c                    1 <= jmu <= ni
c                    jml may be <= 0, and jpl & jpu may be >= ni
c         nphbv = number of elements of arrays jm, jp, km, kp
c         jm = lower point to use
c         jp = upper point to use
c         km = circle at lower limit of segment
c         kp = circle at upper limit of segment
c         phm = azimuthal angle at lower limit of segment
c         php = azimuthal angle at upper limit of segment
c         ph = azimuthal angle at centre of segment
c         dph = azimuthal length of segment between limits
c Return value: -1 = error (tol is too large)
c               0 = segment is not edge of polygon
c               1 = segment is edge of polygon
c               2 = done all segments
c
      if (jpl.ne.0) then
c        done
        if (jpl.eq.jl+ni) goto 220
      endif
c        sin th(i)
      si=sqrt(cmi*(2._10-cmi))
c        dphc = azimuthal angle corresponding to great circle angle tol
      if (tol.gt.PI) goto 300
      dphc=sin(tol/2._10)/si
c        abort if tol/2 exceeds th(i)
      if (dphc.gt.1._10) goto 300
      dphc=2._10*asin(dphc)
c--------first segment: jml <= 1 <= jmu < jpl <= jpu < jml+ni
      if (jpl.eq.0) then
c        lower point: jml to jmu are all at the same azimuth phi
        jml=1
        do j=jml,ni
          jmu=j
          jp(1)=mod(j,ni)+1
          km(1)=(iord(j)+1)/2
          kp(1)=(iord(jp(1))+1)/2
c        lower and upper angles of segment
          phm=phi(1+mod(iord(j)+1,2),km(1))
          php=phi(1+mod(iord(jp(1))+1,2),kp(1))
          dph=php-phm
          if (dph.lt.0._10) dph=dph+TWOPI
          if (dph.gt.dphc) goto 110
        enddo
  110   continue
c        check if lower point goes lower than 1
        do j=ni,jmu+1,-1
          jp(1)=mod(j,ni)+1
          km(1)=(iord(j)+1)/2
          kp(1)=(iord(jp(1))+1)/2
c        lower and upper angles of segment
          phm=phi(1+mod(iord(j)+1,2),km(1))
          php=phi(1+mod(iord(jp(1))+1,2),kp(1))
          dph=php-phm
          if (dph.lt.0._10) dph=dph+TWOPI
          if (dph.gt.dphc) goto 120
          jml=j
        enddo
  120   continue
        if (jml.gt.jmu) jml=jml-ni
c        record range jl to ju of lower point of first segment
        jl=jml
        ju=jmu
c        segment is entire circle
        if (jmu-jml+1.ge.ni) then
c         print *,'gsegij: segment',jml,' to',jmu,
c    *      ' includes all',ni,' intersections'
          jpl=jmu+1
          jpu=jmu+ni
        else
c        upper point: jpl to jpu are all at the same azimuth phi
          jpl=jmu+1
          do j=jpl,jl-1+ni
            jpu=j
            jp(1)=mod(j,ni)+1
            km(1)=(iord(j)+1)/2
            kp(1)=(iord(jp(1))+1)/2
c        lower and upper angles of segment
            phm=phi(1+mod(iord(j)+1,2),km(1))
            php=phi(1+mod(iord(jp(1))+1,2),kp(1))
            dph=php-phm
            if (dph.lt.0._10) dph=dph+TWOPI
            if (dph.gt.dphc) goto 130
          enddo
  130     continue
c        first segment: jml <= 1 <= jmu < jpl <= jpu < jl+ni
          if (jml.gt.1.or.jml.gt.jmu.or.jmu.lt.1.or.jmu.gt.ni
     *      .or.jmu.ge.jpl.or.jpl.gt.jpu.or.jpu.ge.jl+ni) then
c        shouldn't happen
            print *,'*** error from gsegij: 1st segment,',
     *        ni,' intersections:',jml,jmu,jpl,jpu
          endif
        endif
c--------subsequent segments: jml <= jmu < jpl <= jpu <= ju+ni
      elseif (jpl.ne.0) then
c        lower point: jml to jmu are all at the same azimuth phi
        jml=jpl
        do j=jml,jl-1+ni
          jmu=j
          jp(1)=mod(j,ni)+1
          km(1)=(iord(j)+1)/2
          kp(1)=(iord(jp(1))+1)/2
c        lower and upper angles of segment
          phm=phi(1+mod(iord(j)+1,2),km(1))
          php=phi(1+mod(iord(jp(1))+1,2),kp(1))
          dph=php-phm
          if (dph.lt.0._10) dph=dph+TWOPI
          if (dph.gt.dphc) goto 140
        enddo
  140   continue
c        upper point: jpl to jpu are all at the same azimuth phi
        jpl=jmu+1
        do j=jpl,ni+ju
          jpu=j
          jm(1)=mod(j-1+ni,ni)+1
          jp(1)=mod(j,ni)+1
          km(1)=(iord(jm(1))+1)/2
          kp(1)=(iord(jp(1))+1)/2
c        lower and upper angles of segment
          phm=phi(1+mod(iord(jm(1))+1,2),km(1))
          php=phi(1+mod(iord(jp(1))+1,2),kp(1))
          dph=php-phm
          if (dph.lt.0._10) dph=dph+TWOPI
          if (dph.gt.dphc) goto 150
        enddo
  150   continue
c        subsequent segments: 1 < jml <= jmu < jpl <= jpu <= ju+ni
        if (jml.le.1.or.jml.gt.jmu.or.jmu.lt.1.or.jmu.gt.ni
     *    .or.jmu.ge.jpl.or.jpl.gt.jpu.or.jpu.gt.ju+ni) then
c        shouldn't happen
          print *,'*** error from gsegij: >= 2nd segment,',
     *      ni,' intersections:',jml,jmu,jpl,jpu
        endif
      endif
c--------process segment
c        lower angle(s) must be lower limit of segments jml to jmu
      do j=jml,jmu
        jm(1)=mod(j-1+ni,ni)+1
c        lower angle is an upper limit
        if (mod(iord(jm(1)),2).eq.0) then
          km(1)=(iord(jm(1))+1)/2
c        lower angle is not a lower limit
          if (phi(1,km(1)).ne.phi(2,km(1))) goto 210
        endif
      enddo
c        upper angle(s) must be upper limit of segments jpl to jpu
      do j=jpl,jpu
        jp(1)=mod(j-1+ni,ni)+1
c        upper angle is a lower limit
        if (mod(iord(jp(1)),2).eq.1) then
          kp(1)=(iord(jp(1))+1)/2
c        upper angle is not an upper limit
          if (phi(1,kp(1)).ne.phi(2,kp(1))) goto 210
        endif
      enddo
      jm(1)=mod(jml-1+ni,ni)+1
      jp(1)=mod(jpu-1+ni,ni)+1
      km(1)=(iord(jm(1))+1)/2
      kp(1)=(iord(jp(1))+1)/2
      phm=phi(1+mod(iord(jm(1))+1,2),km(1))
      php=phi(1+mod(iord(jp(1))+1,2),kp(1))
c        check segment satisfies all conditions
      if (php.gt.phm) then
        do j=1,ni
c        need to check against k circle only once, at lower limit
          if (mod(iord(j),2).eq.1) then
            k=(iord(j)+1)/2
c        require order   k- ...ph- ph+... k+
            if (phi(2,k).gt.phi(1,k)) then
              if (php.le.phi(1,k).or.phm.ge.phi(2,k)) goto 210
c        or   k+ k- ...ph- ph+...   or   ...ph- ph+... k+ k-
            elseif (phi(2,k).lt.phi(1,k)) then
              if (php.le.phi(1,k).and.phm.ge.phi(2,k)) goto 210
            endif
          endif
        enddo
      elseif (php.lt.phm) then
        do j=1,ni
c        need to check against k circle only once, at lower limit
          if (mod(iord(j),2).eq.1) then
            k=(iord(j)+1)/2
c        require order   ph+... k+ k- ...ph-
            if (phi(2,k).gt.phi(1,k)) then
              if (php.le.phi(1,k).and.phm.ge.phi(2,k)) goto 210
            endif
          endif
        enddo
      endif
c........point at lower limit is one with largest exterior angle psi
      if (jml.eq.jmu) then
        if (nphbv.eq.2) then
          jj=jm(1)
          kk=km(1)
          do iphbv=1,nphbv
            if ((iphbv.eq.1
     *        .and.((i.le.npb.and.kk.le.npb)
     *          .or.(i.gt.npb.and.kk.gt.npb)
     *          .or.kk.gt.npc))
     *      .or.(iphbv.eq.2
     *        .and.((i.le.npb.and.kk.gt.npb)
     *          .or.(i.gt.npb.and.kk.le.npb)
     *          .or.kk.gt.npc))) then
              jm(iphbv)=jj
              km(iphbv)=kk
            else
              jm(iphbv)=0
              km(iphbv)=0
            endif
          enddo
        endif
      else
        do iphbv=1,nphbv
          psim(iphbv)=-1._10-2._10*tol
        enddo
        if (nphbv.eq.2) then
          do iphbv=1,nphbv
            jm(iphbv)=0
            km(iphbv)=0
          enddo
        endif
        do j=jml,jmu
          jj=mod(j-1+ni,ni)+1
c        lower angle is a lower limit
          if (mod(iord(jj),2).eq.1) then
            kk=(iord(jj)+1)/2
            cmk=abs(cm(kk))
c        cmik = 1-cos th(ik)
            cmik=((rpi(1)-rp(1,kk))**2+(rpi(2)-rp(2,kk))**2
     *        +(rpi(3)-rp(3,kk))**2)/2._10
c        bik = cik-ci*ck
c        d = 1-ci^2-ck^2-cik^2+2*ci*ck*cik
c        cos psi = bik/(si*sk)
c        sin psi = sqrt(d)/(si*sk)
c        psi = atan(sqrt(d)/bik) is exterior angle at intersection
            bik=(cmi+cmk)-cmi*cmk-cmik
            d=-(cmi-cmk)**2+cmik*(2._10*((cmi+cmk)-cmi*cmk)-cmik)
            if (d.lt.0._10) d=0._10
            if ((scmi.ge.0.and.cm(kk).lt.0._10)
     *        .or.(scmi.lt.0.and.cm(kk).ge.0._10)) bik=-bik
            do iphbv=1,nphbv
              if (iphbv.eq.2) bik=-bik
              psi=atan2(sqrt(d),bik)
              ismax=.false.
              if (nphbv.eq.1
     *          .or.(iphbv.eq.1
     *            .and.((i.le.npb.and.kk.le.npb)
     *              .or.(i.gt.npb.and.kk.gt.npb)
     *              .or.kk.gt.npc))
     *          .or.(iphbv.eq.2
     *            .and.((i.le.npb.and.kk.gt.npb)
     *              .or.(i.gt.npb.and.kk.le.npb)
     *              .or.kk.gt.npc))) then
c        choose largest exterior angle psi
                if (psi.gt.psim(iphbv)+tol) then
                  ismax=.true.
                elseif (psi.ge.psim(iphbv)-tol) then
                  if (cm(kk).ge.0._10) then
                    if (cm(km(iphbv)).ge.0._10) then
c        only do tighter of two circles with same exterior angle
                      if (cm(kk).lt.cm(km(iphbv))) then
                        ismax=.true.
c        or later of two degenerate circles
                      elseif (cm(kk).eq.cm(km(iphbv))
     *                  .and.kk.gt.km(iphbv)) then
                        ismax=.true.
                      endif
                    else
                      if (cm(kk)-1._10.lt.cm(km(iphbv))+1._10) then
                        ismax=.true.
                      elseif (cm(kk)-1._10.eq.cm(km(iphbv))+1._10
     *                  .and.kk.gt.km(iphbv)) then
                        ismax=.true.
                      endif
                    endif
                  else
                    if (cm(km(iphbv)).lt.0._10) then
                      if (cm(kk).lt.cm(km(iphbv))) then
                        ismax=.true.
                      elseif (cm(kk).eq.cm(km(iphbv))
     *                  .and.kk.gt.km(iphbv)) then
                        ismax=.true.
                      endif
                    else
                      if (cm(kk)+1._10.lt.cm(km(iphbv))-1._10) then
                        ismax=.true.
                      elseif (cm(kk)+1._10.eq.cm(km(iphbv))-1._10
     *                  .and.kk.gt.km(iphbv)) then
                        ismax=.true.
                      endif
                    endif
                  endif
                endif
              endif
              if (ismax) then
                jm(iphbv)=jj
                km(iphbv)=kk
                psim(iphbv)=psi
              endif
            enddo
          endif
        enddo
      endif
c........point at upper limit is one with largest exterior angle psi
      if (jpl.eq.jpu) then
        if (nphbv.eq.2) then
          jj=jp(1)
          kk=kp(1)
          do iphbv=1,nphbv
            if ((iphbv.eq.1
     *        .and.((i.le.npb.and.kk.le.npb)
     *          .or.(i.gt.npb.and.kk.gt.npb)
     *          .or.kk.gt.npc))
     *      .or.(iphbv.eq.2
     *        .and.((i.le.npb.and.kk.gt.npb)
     *          .or.(i.gt.npb.and.kk.le.npb)
     *          .or.kk.gt.npc))) then
              jp(iphbv)=jj
              kp(iphbv)=kk
            else
              jp(iphbv)=0
              kp(iphbv)=0
            endif
          enddo
        endif
      else
        do iphbv=1,nphbv
          psim(iphbv)=-1._10-2._10*tol
        enddo
        if (nphbv.eq.2) then
          do iphbv=1,nphbv
            jp(iphbv)=0
            kp(iphbv)=0
          enddo
        endif
        do j=jpl,jpu
          jj=mod(j-1+ni,ni)+1
c        upper angle is an upper limit
          if (mod(iord(jj),2).eq.0) then
            kk=(iord(jj)+1)/2
            cmk=abs(cm(kk))
c        cmik = 1-cos th(ik)
            cmik=((rpi(1)-rp(1,kk))**2+(rpi(2)-rp(2,kk))**2
     *        +(rpi(3)-rp(3,kk))**2)/2._10
c        bik = cik-ci*ck
c        d = 1-ci^2-ck^2-cik^2+2*ci*ck*cik
c        cos psi = bik/(si*sk)
c        sin psi = sqrt(d)/(si*sk)
c        psi = atan(sqrt(d)/bik) is exterior angle at intersection
            bik=(cmi+cmk)-cmi*cmk-cmik
            d=-(cmi-cmk)**2+cmik*(2._10*((cmi+cmk)-cmi*cmk)-cmik)
            if (d.lt.0._10) d=0._10
            if ((scmi.ge.0.and.cm(kk).lt.0._10)
     *        .or.(scmi.lt.0.and.cm(kk).ge.0._10)) bik=-bik
            do iphbv=1,nphbv
              if (iphbv.eq.2) bik=-bik
              psi=atan2(sqrt(d),bik)
              ismax=.false.
              if (nphbv.eq.1
     *          .or.(iphbv.eq.1
     *            .and.((i.le.npb.and.kk.le.npb)
     *              .or.(i.gt.npb.and.kk.gt.npb)
     *              .or.kk.gt.npc))
     *          .or.(iphbv.eq.2
     *            .and.((i.le.npb.and.kk.gt.npb)
     *              .or.(i.gt.npb.and.kk.le.npb)
     *              .or.kk.gt.npc))) then
c        choose largest exterior angle psi
                if (psi.gt.psim(iphbv)+tol) then
                  ismax=.true.
                elseif (psi.ge.psim(iphbv)-tol) then
                  if (cm(kk).ge.0._10) then
                    if (cm(kp(iphbv)).ge.0._10) then
c        only do tighter of two circles with same exterior angle
                      if (cm(kk).lt.cm(kp(iphbv))) then
                        ismax=.true.
c        or later of two degenerate circles
                      elseif (cm(kk).eq.cm(kp(iphbv))
     *                  .and.kk.gt.kp(iphbv)) then
                        ismax=.true.
                      endif
                    else
                      if (cm(kk)-1._10.lt.cm(kp(iphbv))+1._10) then
                        ismax=.true.
                      elseif (cm(kk)-1._10.eq.cm(kp(iphbv))+1._10
     *                  .and.kk.gt.kp(iphbv)) then
                        ismax=.true.
                      endif
                    endif
                  else
                    if (cm(kp(iphbv)).lt.0._10) then
                      if (cm(kk).lt.cm(kp(iphbv))) then
                        ismax=.true.
                      elseif (cm(kk).eq.cm(kp(iphbv))
     *                  .and.kk.gt.kp(iphbv)) then
                        ismax=.true.
                      endif
                    else
                      if (cm(kk)+1._10.lt.cm(kp(iphbv))-1._10) then
                        ismax=.true.
                      elseif (cm(kk)+1._10.eq.cm(kp(iphbv))-1._10
     *                  .and.kk.gt.kp(iphbv)) then
                        ismax=.true.
                      endif
                    endif
                  endif
                endif
              endif
              if (ismax) then
                jp(iphbv)=jj
                kp(iphbv)=kk
                psim(iphbv)=psi
              endif
            enddo
          endif
        enddo
      endif
c        circles at upper and lower limits
      if (km(1).ne.0) then
        phm=phi(1+mod(iord(jm(1))+1,2),km(1))
      else
        phm=phi(1+mod(iord(jm(2))+1,2),km(2))
      endif
      if (kp(1).ne.0) then
        php=phi(1+mod(iord(jp(1))+1,2),kp(1))
      else
        php=phi(1+mod(iord(jp(2))+1,2),kp(2))
      endif
c        angular length, centre point of segment
      if (php.gt.phm) then
        dph=php-phm
        ph=(php+phm)/2._10
      elseif (php.le.phm) then
        ph=(php+phm)/2._10
        if (ph.le.0._10) then
          ph=ph+PI
          php=php+TWOPI
        elseif (ph.gt.0._10) then
          ph=ph-PI
          phm=phm-TWOPI
        endif
        dph=php-phm
      endif
c      segment is edge of polygon
  200 gsegij=1
      return
c
c------segment is not edge of polygon
  210 gsegij=0
      return
c
c------done full circle
  220 gsegij=2
      return
c
c------tol is too large
  300 gsegij=-1
      return
c
      end
c
c-----------------------------------------------------------------------
      subroutine gvtrail(scmi,np,i,km,kp,vtrail,vik,nvmax,nv,ik,ikchk)
      integer scmi,np,i,km,kp,nvmax,vtrail(nvmax),vik(nvmax,2),nv,ik
      real*10 ikchk
c
c        externals
      integer ik2ik
c     integer ik2i,ik2k
c        local (automatic) variables
      integer iv,k,l
      real*10 ikran
c *
c * Record endpoints of edge.
c *
c  Input: scmi
c         np
c         i
c         km, kp
c         nvmax
c         nv
c Output: vtrail
c         vik
c         ik
c Input/Output: ikchk
c         
c        two end points of the edge, going right-handedly around edge
      do l=1,2
c        end point is intersection of i circle with k circle
        if ((l.eq.1.and.scmi.ge.0).or.(l.eq.2.and.scmi.lt.0)) then
          k=km
        else
          k=kp
        endif
c       print *,'     intersect',i,k
c        from edge k to edge i right-handedly through vertex
        if (l.eq.1) then
          ik=ik2ik(np,k,i)
c        from edge i to edge k right-handedly through vertex
        else
          ik=ik2ik(np,i,k)
        endif
        if (nv.le.nvmax) then
          vik(nv,l)=ik
        endif
c        pseudo-random number from ik
        call ikrand(ik,ikran)
c        new intersection
        if (i.lt.k) then
c        ikchk = ikchk + ikran, added as unsigned long long's
          call ikrandp(ikchk,ikran)
c        intersection already met
        else
c        ikchk = ikchk - ikran, subtracted as unsigned long long's
          call ikrandm(ikchk,ikran)
          if (nv.le.nvmax) then
c        index of vertex already met
            if (l.eq.1) then
              do iv=1,nv
                if (vik(iv,2).eq.ik) then
c        vertex after iv is nv
                  vtrail(iv)=nv
                  goto 240
                endif
              enddo
            elseif (l.eq.2) then
              do iv=1,nv
                if (vik(iv,1).eq.ik) then
c        vertex after nv is iv
                  vtrail(nv)=iv
                  goto 240
                endif
              enddo
            endif
c        can happen if near multiple intersection is inconsistent
c           print *,'*** from gvtrail: no vertex',
c    *        ' (',ik2i(np,ik),ik2k(np,ik),')',
c    *        ' found at',i,', edge',km,kp,'???'
c           write (*,'(" vertices so far are:",$)')
c           do iv=1,nv
c             write (*,'(" (",i2,i3,") (",i2,i3,")",$)')
c    *          ik2k(np,vik(iv,1)),ik2i(np,vik(iv,1)),
c    *          ik2i(np,vik(iv,2)),ik2k(np,vik(iv,2))
c           enddo
c           write (*,'(/$)')
  240       continue
          endif
        endif
      enddo
      return
      end
c
c-----------------------------------------------------------------------
      integer function ik2ik(np,i,k)
      integer np,i,k
c
      integer ik
c
c        gives ik = 1,2,...,np(np-1) for i,k = 0,1,...,np with i != k
      ik=i+np*k
      if (i.lt.k) ik=ik+1
      ik2ik=ik
      return
      end
c
c-----------------------------------------------------------------------
      integer function ik2i(np,ik)
      integer np,ik
c
      integer i,k
c
      i=mod(ik-1,np)
      k=(ik-1)/np
      if (i.ge.k) i=i+1
      ik2i=i
      return
      end
c
c-----------------------------------------------------------------------
      integer function ik2k(np,ik)
      integer np,ik
c
      integer k
c
      k=(ik-1)/np
      ik2k=k
      return
      end
c
c-----------------------------------------------------------------------
      subroutine gvord(np,nv,ipv,gp,ev,nev,vtrail,vord,vik,gord)
      integer np,nv,ipv(nv),gp(np),ev(nv),nev,
     *  vtrail(nv),vord(nv),vik(nv),gord(nv)
c
c        externals
      integer ik2i,ik2k
c        local (automatic) variables
      integer iev,inter,iv,ivm,jv,kv,nvm
c *
c * Order of vertices around polygon.
c *
c  Input: np = number of caps
c         nv = number of vertices
c         ipv(iv) = circle number of vertex iv
c         gp(i) = group to which circle i belongs
c         vtrail = vertex trail
c         vik(iv) = value of ik at iv'th vertex
c Output: vord = vertex order
c         ev = end indices of each connected sequence of vertices
c         nev = number of connected sequences of vertices
c
      nev=0
      if (nv.gt.0) then
c       print *,'gvord: vertex trail:',(vtrail(iv),iv=1,nv)
c       print *,'gvord: vik_1 on input:',(vord(iv),iv=1,nv)
c       print *,'gvord: vik_2 on input:',(vik(iv),iv=1,nv)
c        group vertex belongs to
        do iv=1,nv
          vord(iv)=gp(ipv(iv))
        enddo
c        order vertices by group
        call finibot(vord,nv,gord,nv)
c        initialize vord to zero
        do iv=1,nv
          vord(iv)=0
        enddo
        jv=0
        iev=0
c        do intersecting, then non-intersecting circles
        do inter=1,2
          do 320 ivm=1,nv
c        order connected sequences by group
            nvm=gord(ivm)
c        find first vertex not already traversed
            do iv=1,nv
              if (vord(iv).eq.nvm) goto 320
            enddo
            if (inter.eq.1) then
              if (ik2i(np,vik(nvm)).eq.0.or.ik2k(np,vik(nvm)).eq.0)
     *          goto 320
            else
              if (ik2i(np,vik(nvm)).ne.0.and.ik2k(np,vik(nvm)).ne.0)
     *          goto 320
            endif
c        index of first vertex of this boundary
            jv=jv+1
            vord(jv)=nvm
c        circulate
            kv=jv+1
            do jv=kv,nv
c        next vertex of this boundary
              vord(jv)=vtrail(vord(jv-1))
c        gone full circle
              if (vtrail(vord(jv)).eq.nvm) then
c        record cumulative length
                iev=iev+1
                ev(iev)=jv
                nev=nev+1
                if (jv.eq.nv) goto 330
                goto 320
              endif
            enddo
  320     continue
        enddo
  330   continue
c       print *,'vertex order:',(vord(iv),iv=1,nv)
c       print *,'vertex groups:',(gp(ipv(iv)),iv=1,nv)
c        should not happen
        if (vtrail(vord(nv)).ne.nvm) then
          print *,'*** from gvord: last vertex',vord(nv),
     *      ' should connect to',nvm,'???'
        endif
      endif
      return
c
      end
c
c-----------------------------------------------------------------------
      real*10 function cmijf(rpi,rpj)
      real*10 rpi(3),rpj(3)
c *
c * 1 - cos th(ij)
c * where th(ij) is angle between unit vectors rpi and rpj.
c *
c * Input: rpi, rpj = unit vectors
c  Output: cmij = 1 - cos th(ij)
c
      cmijf=((rpi(1)-rpj(1))**2+(rpi(2)-rpj(2))**2
     *  +(rpi(3)-rpj(3))**2)/2._10
c
      return
      end
c
c-----------------------------------------------------------------------
      subroutine vpermi(ix,n,iperm,iwk)
      integer n,iperm(n)
      integer ix(n),iwk(n)
c
c        local (automatic) variables
      integer i,j
c *
c * Permutes elements of integer array ix(n) by permutation iperm(n).
c *
c  Input: n = dimension of array ix
c         iperm = array of dimension n specifying permutation of ix
c Input/Output: ix = integer array of dimension n
c Work array: wk of dimension n
c
      do i=1,n
        iwk(i)=ix(i)
      enddo
      do i=1,n
        j=iperm(i)
        ix(i)=iwk(j)
      enddo
      return
      end
c
c-----------------------------------------------------------------------
      subroutine vpermd(x,n,iperm,wk)
      integer n,iperm(n)
      real*10 x(n),wk(n)
c
c        local (automatic) variables
      integer i,j
c *
c * Permutes elements of real*10 array x(n) by permutation iperm(n).
c *
c  Input: n = dimension of array x
c         iperm = array of dimension n specifying permutation of x
c Input/Output: x = real*10 array of dimension n
c Work array: wk of dimension n
c
      do i=1,n
        wk(i)=x(i)
      enddo
      do i=1,n
        j=iperm(i)
        x(i)=wk(j)
      enddo
      return
      end
c
c-----------------------------------------------------------------------
      subroutine vpermdd(x,m,n,iperm,wk)
      integer m,n,iperm(n)
      real*10 x(m,n),wk(n)
c
c        local (automatic) variables
      integer i,j,k
c *
c * Permutes each column of real*10 array x(m,n) by permutation iperm(n).
c *
c  Input: m,n = dimensions of array x
c         iperm = array of dimension n specifying permutation
c                 of each column x
c Input/Output: x = real*10 array of dimension m,n
c Work array: wk of dimension n
c
      do k=1,m
        do i=1,n
          wk(i)=x(k,i)
        enddo
        do i=1,n
          j=iperm(i)
          x(k,i)=wk(j)
        enddo
      enddo
      return
      end
c
c-----------------------------------------------------------------------
      integer function garpi(area,iarea,rp,cm,np,
     *  whole,nbd0m,nbd0p,nbd,nmult,tol)
      integer iarea,np,nbd0m,nbd0p,nbd,nmult
      logical whole
      real*10 area,rp(3,np),cm(np),tol
c
c        parameters
      include 'pi.par'
      real*10 TWOPI
      parameter (TWOPI=2._10*PI)
c        intrinsics
      intrinsic abs
c        data variables
      real*10 areatol
c        local (automatic) variables
      integer i,icmmin
      real*10 cmmin,darea,p
c *
c * Add iarea*2*pi to area.
c *
c  Input: rp
c         cm
c         np
c         whole = whether region is whole sphere.
c         nbd0m = number of non-intersecting circles bounding polygon
c                 with cm < 0.
c         nbd0p = number of non-intersecting circles bounding polygon
c                 with cm >= 0.
c         nbd = number of edges bounding polygon,
c               excluding non-intersecting circles.
c         nmult = number of near multiply-intersecting circles
c                 on boundary of polygon.
c Output: iarea = number of 2*pi's by which area was adjusted.
c Input/Output: area -> area+iarea*2*pi.
c Return value: 0 = ok;
c               1 = recommend retry with enlarged tol.
c
c        ok if area tests not too far outside [0,max]
      data areatol /1.e-10_10/
c
      if (whole) then
        iarea=2
        area=area+iarea*TWOPI
c        all boundaries are non-intersecting circles
      elseif (nbd.eq.0) then
        if (nbd0m.eq.0.and.nbd0p.eq.0) then
          iarea=0
        else
c        from area formula involving the Euler characteristic
          iarea=2*(1-nbd0p)
          area=area+iarea*TWOPI
        endif
c        some boundaries are intersecting circles
      else
c        add/subtract 2*pi's to area to ensure 0.le.area.lt.2*pi
        iarea=area/TWOPI
        area=area-iarea*TWOPI
C       if (iarea.ne.0) print *,'area +=',iarea,' * TWOPI =',area
        if (area.lt.0._10) then
          iarea=iarea+1
          area=area+TWOPI
C         print *,'area += TWOPI =',area
        endif
c        there were near multiple intersections
        if (nmult.ge.1) then
c        chances are area just less than 2*pi is actually zero
          if (area.ge.TWOPI-areatol) then
            iarea=iarea-1
            area=0._10
            goto 400
          endif
        endif
c        check area does not exceed area within any one circle
        cmmin=2._10
        do i=1,np
          if (cm(i).ge.0._10) then
            if (cm(i).lt.cmmin) then
              cmmin=cm(i)
              icmmin=i
            endif
          elseif (cm(i).lt.0._10) then
            if (2._10+cm(i).lt.cmmin) then
              cmmin=2._10+cm(i)
              icmmin=i
            endif
          endif
        enddo
        darea=TWOPI*cmmin
        p=area/darea
C       print *,'area/area(',icmmin,')=',area,' /',darea,' =',p
        if (p.gt.1._10) then
c        check if discrepancy is from numerical roundoff
          if (abs(area-TWOPI).le.areatol) then
            area=0._10
          elseif (area.le.darea+areatol) then
            area=darea
c        problem is genuine: can happen with nearly kissing circles
          else
C           print *,'*** from garpi: at tol =',tol,
C    *        ' area/area(',icmmin,')=',area,' /',darea,
C    *        ' =',p,'  should be .le. 1'
            goto 410
          endif
        endif
      endif
c        done
 400  garpi=0
      return
c
c        return with recommendation to retry with enlarged tol
 410  garpi=1
      return
      end
c
c-----------------------------------------------------------------------
      subroutine gtol(tol,tolin)
      real*10 tol,tolin
c *
c * Modify tolerance tol to multiple intersections.
c * The tolerance tol is changed by successive factors of 2
c * from the original input tolerance tolin,
c * see-sawing between smaller and larger values.
c *
      if (tolin.le.0._10) then
c        write(*,*) 'tolin =', tolin
        if (tol.le.0._10) then
          tol=1.e-15_10
c          write(*,*) 'tol =', tol
        else
          tol=tol*4._10
c          write(*,*) 'tol =', tol
        endif
c        see-saw tolerance between smaller and larger values
      else
        if (tol.ge.tolin) then
c          write(*,*) '(GE) before: TOL=', tol
c          write(*,*) '(ge) before: tolin=', tolin
          tol=tolin*tolin/tol/6._10
c          if (tol.gt.0.00001) then
c            write(*,*) '(ge) after: tol=', tol
c          endif
c          write(*,*) '(ge) after: tolin=', tolin
        else
c          write(*,*) '(lt) before: tol=', tol
c          write(*,*) '(lt) before: tolin=', tolin
          tol=tolin*tolin/tol
c          if (tol.ge.1.e-5_10) then
c            write(*,*) '(lt) after: tol=', tol
c          endif
c          write(*,*) '(lt) after: tolin=', tolin
        endif
      endif
      return
      end
c
