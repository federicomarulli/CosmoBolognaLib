c-----------------------------------------------------------------------
c © A J S Hamilton 2001
c-----------------------------------------------------------------------
c * This routine needs upgrading to the same level of robustness
c * as the rest of the mangle software:
c * 1. Should check consistency of multiple intersections,
c *    and increase tol and repeat if an inconsistency is detected.
c * 2. As a corollary of 1., need to go around all edges,
c *    not just the particular edge i.
c * 3. If 2. is done, then the routine should be able to deal
c *    with point abuts in addition to abutments of finite extent.
c *
c * Currently, the routine deals correctly with multiple intersections,
c * albeit without checking for consistency.
c * However, it deals only with abutments of finite extent,
c * NOT with point abutments, where polygons abut at a single point
c * or at a set of isolated points.
c * This is enough to get bound right,
c * but vert is missing contributions from point abuts.
c-----------------------------------------------------------------------
      subroutine gphbv(bound,vert,rp,cm,np,npb,npc,i,tol,phi,iord)
      integer np,npb,npc,i,iord(2*np)
      real*10 bound(2),vert(2),rp(3,np),cm(np),tol,phi(2,np)
c
c        parameters
      include 'pi.par'
      real*10 TWOPI
      parameter (TWOPI=2._10*PI)
c        intrinsics
      intrinsic abs
c        externals
      integer gsegij,gzeroar
c        data variables
      real*10 big,bndtol,psitol
      real*10 dphmin
c        local (automatic) variables
      integer iphbv,iseg,j,jm(2),jml,jmu,jp(2),jpl,jpu,k,km(2),kp(2),
     *  l,ni,nmult
      logical warn
      real*10 bik,cmi,cmik,cmk,cti(3),ctk(3),ctpsi(3),
     *  d,dbound(2),dph,dvert(2),p,ph,phm,php,psi(3),
     *  scmi,si,sk,t(3),xi(3),yi(3)
c *
c * Correction to boundary and vertex terms per subroutine gspher
c * arising from disjoint regions W12 and W13 which abut along circle i
c * defined by 1 - r.rp(i) = cm(i).
c * W12 (W13) is the intersection of W2 (W3) with global region W1.
c * Regions are bounded by
c *    1 - r.rp(j) < cm(j)  (if cm(j).ge.0)
c *    1 - r.rp(j) > -cm(j)  (if cm(j).lt.0)
c * for j=1,np where rp(j) are unit directions.
c * The boundary of W2 is in j =     1 to npb
c *                 W3           npb+1 to npc
c *                 W1           npc+1 to np.
c * If circle i belongs to W2 (i.e. i=1 to npb) then
c * W2 (W3) is north (south) of circle i if cm(i).ge.0,
c *            south (north) of circle i if cm(i).lt.0.
c * Conversely, if circle i belongs to W3 (i.e. i=npb+1 to npc) then
c * W3 (W2) is north (south) of circle i if cm(i).ge.0,
c *            south (north) of circle i if cm(i).lt.0.
c * The abutting boundary i' of W3 (W2) which coincides with
c * boundary i of W2 (W3) should be suppressed by setting cm(i')=2.
c * Boundaries j that are duplicates of boundary i
c * or j' of the abutting boundary i'
c * should be suppressed by setting cm(j)=2 and cm(j')=2.
c *
c * Boundary and vertex terms are evaluated for -<W12 W13>.
c *
c * To get bound & vert for correlation <(W1-W12-W13)^2>
c * of region W1 less its intersection with disjoint regions W2 and W3,
c * follow instructions in gspher, i.e. basically
c *    call gspher(area,bound,vert,parameters of W1,ibv=0,...)
c *    call gspher(darea,dbound,dvert,parameters of W2 & W1,ibv=2,...)
c *    area=area-darea, etc.
c *    call gspher(darea,dbound,dvert,parameters of W3 & W1,ibv=2,...)
c *    area=area-darea, etc.
c * and then, if W2 and W3 abut,
c *    call gphbv(dbound,dvert,parameters of W2 W3 & W1,...)
c *    bound=bound-2*dbound
c *    vert=vert-2*dvert
c *
c * To get bound & vert for cross-correlation <W12(W1-W13)>
c * between W12
c * and region W1 less its intersection with W3 outside W2,
c * follow instructions in gspher, i.e.
c *    call gspher(area,bound,vert,parameters of W2 & W1,ibv=1,...)
c * and then, if W2 and W3 abut,
c *    call gphbv(dbound,dvert,parameters of W2 W3 & W1,...)
c *    bound=bound+dbound
c *    vert=vert+dvert
c *
c * To get bound & vert for cross-correlation <(W12-W13)(W1-W13)>
c * between W12 less its intersection with W3 inside W2,
c * and region W1 less its intersection with W3,
c * follow instructions in gspher, i.e.
c *    call gspher(area,bound,vert,parameters of W2 & W1,ibv=1,...)
c *    call gspher(darea,dbound,dvert,parameters of W3 & W12,ibv=2,...)
c *    area=area-darea, etc.
c * and then, if W2' and W3 abut, where W2' denotes complement of W2,
c *    call gphbv(dbound,dvert,parameters of W2' W3 & W1,...)
c *    bound=bound-dbound
c *    vert=vert-dvert
c * Note that usually all W2' constraints should be null.
c * If however W2' abuts W3 along two distinct boundaries, then W2'
c * should be split into two regions W2'a and W2'b which each have a
c * single abutment with W3, and two calls to gphbv should be made
c *    call gphbv(dbound,dvert,parameters of W2'a W3 & W1,...)
c *    bound=bound-dbound
c *    vert=vert-dvert
c *    call gphbv(dbound,dvert,parameters of W2'b W3 & W1,...)
c *    bound=bound-dbound, etc.
c *
c * In general, to get bound & vert for cross-correlation
c * <(W12-W13-...)(W1-W13-W14-...)>
c * between region W12 less its intersection W13 with a
c * bunch of disjoint W3 lying inside W2,
c * and W1 less its intersections W13 and W14 with that bunch of W3 and
c * another bunch of disjoint W4 lying outside W2,
c *    call gspher(area,bound,vert,parameters of W2 & W1,ibv=1,...)
c * then chop out all the W13s by calls
c *    call gspher(darea,dbound,dvert,parameters of W3 & W12,ibv=2,...)
c *    area=area-darea, etc.
c * then for each W3 (inside W2) abutting the edge of W2
c *    call gphbv(dbound,dvert,parameters of W2' W3 & W1,...)
c *    bound=bound-dbound, etc.
c * for each W4 (outside W2) abutting W2
c *    call gphbv(dbound,dvert,parameters of W2 W4 & W1,...)
c *    bound=bound+dbound, etc.
c * for each pair W3a & W3b (inside W2) which abut each other
c *    call gphbv(dbound,dvert,parameters of W3a W3b & W1,...)
c *    bound=bound-2*dbound, etc.
c * for each pair W3 (inside W2) & W4 (outside W2) abutting each other
c *    call gphbv(dbound,dvert,parameters of W3 W4 & W1,...)
c *    bound=bound-dbound, etc.
c *
c  Input: rp(3,j),j=1,np
c         cm(j),j=1,np
c         np, npb, npc: W2     1 to npb
c                       W3 npb+1 to npc
c                       W1 npc+1 to np.
c         i = abutting boundary of W2 (or W3);
c             the abutting boundary i' of W3 (or W2) which coincides
c             with boundary i of W2 (or W3) should be suppressed
c             by setting cm(i')=2.
c Output: bound(2)
c         vert(2)
c Work arrays: phi and iord should be dimensioned at least 2*np
c
c        set azimuthal angle of non-intersection to big
      data big /1.e6_10/
c        set vertex term to zero if |psi| < psitol
      data psitol /1.e-10_10/
c        ok if bound(1) tests not too far outside [0,max]
      data bndtol /1.e-10_10/
c        warn about multiple intersection when dph < dphmin
      data dphmin /1.e-8_10/
c
C     print *,'--------------------'
c        abutting boundary must belong to W2 or W3
      if (i.gt.npc) then
        print *,'*** from gphbv: i =',i,' should be .le.',npc
        goto 410
      endif
c        zero stuff
      bound(1)=0._10
      bound(2)=0._10
      vert(1)=0._10
      vert(2)=0._10
      warn=.false.
c        check for zero angle because one circle is null
      if (gzeroar(cm,np).eq.0) goto 410
c        cm(i).ge.2 means include whole sphere, which is no constraint
      if (cm(i).ge.2._10) goto 410
c--------identify boundary segments around circle i
      if (cm(i).ge.0._10) then
        scmi=1
      else
        scmi=-1
      endif
      cmi=abs(cm(i))
      si=sqrt(cmi*(2._10-cmi))
c........construct cartesian axes with z-axis along rp(i)
      call gaxisi(rp(1,i),xi,yi)
c........angles phi about z-axis rp(i) of intersection of i & j circles
      call gphij(rp,cm,np,i,rp(1,i),scmi,cmi,xi,yi,big,tol,ni,phi)
c        i circle lies outside polygon
      if (ni.eq.-1) goto 410
c        area of polygon is zero
      if (ni.eq.-2) goto 410
c........i circle has no intersections
      if (ni.eq.0) then
        dph=TWOPI
        dbound(1)=si*dph
        dbound(2)=(1._10/si-2._10*si)*dph
        bound(1)=bound(1)+dbound(1)
        bound(2)=bound(2)+dbound(2)
C       print *,'full circle'
C       print *,'dbound =',dbound(1),dbound(2),
C    *    ' bound =',bound(1),bound(2)
c........i circle has intersections
      elseif (ni.gt.0) then
c        find ordering of intersection angles around i circle
        call findbot(phi,2*np,iord,ni)
c........contribution from each segment of i circle
        jpl=0
c        come here to do another segment
  220   continue
c........is segment edge of polygon?
          iseg=gsegij(rp,cm,np,npb,npc,i,rp(1,i),scmi,cmi,tol,ni,
     *      phi,iord,jml,jmu,jpl,jpu,2,jm,jp,km,kp,phm,php,ph,dph)
c        error
          if (iseg.eq.-1) goto 420
c        not an edge
          if (iseg.eq.0) goto 220
c        gone full circle
          if (iseg.eq.2) goto 280
c        near multiple intersection
          if (dph.lt.dphmin) then
c        increment count of near multiple intersections
            nmult=nmult+1
c           warn=.true.
c           print *,
c    *        '*** warning from gspher: near multiple intersection at'
c    *        ,i,': segment',km(1),kp(1),' &',km(2),kp(2),' dph=',dph
          endif
c........segment satisfies conditions
c. . . . boundary terms
          dbound(1)=si*dph
          dbound(2)=(1._10/si-2._10*si)*dph
          bound(1)=bound(1)+dbound(1)
          bound(2)=bound(2)+dbound(2)
C         print *,'at',i,': edge',km(1),kp(1),' &',km(2),kp(2),
C    *     ' (',jm(1),' &',jm(2),' in',jml,jmu,',',
C    *     jp(1),' &',jp(2),' in',jpl,jpu,' of',ni,')'
C         print *,'dph/(2*pi)=',dph/TWOPI
C         print *,'dbound =',dbound(1),dbound(2),
C    *      ' bound =',bound(1),bound(2)
c. . . . vertex terms
          do l=1,2
            do iphbv=1,2
c        end point is intersection of i circle with k circle
              if (l.eq.1) then
                k=km(iphbv)
              elseif (l.eq.2) then
                k=kp(iphbv)
              endif
              if (k.eq.0) then
                psi(iphbv)=0._10
                ctpsi(iphbv)=1._10/psi(iphbv)
                t(iphbv)=0._10
c        cti = cot th(i)
                cti(iphbv)=(1._10-cmi)/si
                if (scmi.lt.0) cti(iphbv)=-cti(iphbv)
                if (iphbv.eq.2) cti(iphbv)=-cti(iphbv)
c        ctk = cot th(k)
                ctk(iphbv)=cti(iphbv)
              else
                cmk=abs(cm(k))
                sk=sqrt(cmk*(2._10-cmk))
c        cmik = 1-cos th(ik)
                cmik=((rp(1,i)-rp(1,k))**2+(rp(2,i)-rp(2,k))**2
     *            +(rp(3,i)-rp(3,k))**2)/2._10
c        bik = cik-ci*ck
c        d = 1-ci^2-ck^2-cik^2+2*ci*ck*cik
c        cos psi = bik/(si*sk)
c        sin psi = sqrt(d)/(si*sk)
c        psi = atan(sqrt(d)/bik) is exterior angle at intersection
                bik=(cmi+cmk)-cmi*cmk-cmik
                if ((scmi.ge.0.and.cm(k).lt.0._10)
     *            .or.(scmi.le.0.and.cm(k).ge.0._10)) bik=-bik
                if (iphbv.eq.2) bik=-bik
c        i and k circles kiss
                if (phi(1,k).eq.phi(2,k)) then
                  d=0._10
                else
                  d=-(cmi-cmk)**2+cmik*(2._10*((cmi+cmk)-cmi*cmk)-cmik)
c        assert that circles at least touch
                  if (d.lt.0._10) d=0._10
                  d=sqrt(d)
                endif
                ctpsi(iphbv)=bik/d
                psi(iphbv)=atan2(d,bik)
c        t=tan psi/2
                if (bik.gt.0._10) then
                  t(iphbv)=d/(bik+sqrt(bik**2+d**2))
                elseif (bik.lt.0._10) then
                  t(iphbv)=(-bik+sqrt(bik**2+d**2))/d
                elseif (bik.eq.0._10) then
                  t(iphbv)=1._10
                endif
c        cti = cot th(i)
                cti(iphbv)=(1._10-cmi)/si
                if (scmi.lt.0) cti(iphbv)=-cti(iphbv)
                if (iphbv.eq.2) cti(iphbv)=-cti(iphbv)
c        ctk = cot th(k)
                ctk(iphbv)=(1._10-cmk)/sk
                if (cm(k).lt.0._10) ctk(iphbv)=-ctk(iphbv)
              endif
            enddo
c        psi(3) = psi(1) + psi(2) - pi
            psi(3)=psi(1)+psi(2)-PI
c        cot psi(3)
            if (abs(ctpsi(1)).le.1._10) then
              if (abs(ctpsi(2)).le.1._10) then
                ctpsi(3)=(ctpsi(1)*ctpsi(2)-1._10)
     *            /(ctpsi(1)+ctpsi(2))
              else
                ctpsi(3)=(ctpsi(1)-1._10/ctpsi(2))
     *            /(ctpsi(1)/ctpsi(2)+1._10)
              endif
            else
              if (abs(ctpsi(2)).le.1._10) then
                ctpsi(3)=(ctpsi(2)-1._10/ctpsi(1))
     *            /(1._10+ctpsi(2)/ctpsi(1))
              else
                ctpsi(3)=(1._10-1._10/ctpsi(1)/ctpsi(2))
     *            /(1._10/ctpsi(2)+1._10/ctpsi(1))
              endif
            endif
c        tan psi(3)/2
            if (abs(t(1)).le.1._10) then
              if (abs(t(2)).le.1._10) then
                t(3)=(t(1)*t(2)-1._10)/(t(1)+t(2))
              else
                t(3)=(t(1)-1._10/t(2))/(t(1)/t(2)+1._10)
              endif
            else
              if (abs(t(2)).le.1._10) then
                t(3)=(t(2)-1._10/t(1))/(1._10+t(2)/t(1))
              else
                t(3)=(1._10-1._10/t(1)/t(2))/(1._10/t(2)+1._10/t(1))
              endif
            endif
            cti(3)=ctk(1)
            ctk(3)=ctk(2)
            do iphbv=1,3
              if (abs(psi(iphbv)).le.psitol) then
                dvert(1)=0._10
                dvert(2)=0._10
              else
                dvert(1)=1._10-psi(iphbv)*ctpsi(iphbv)
                dvert(2)=t(iphbv)*(3._10+t(iphbv)**2)
     *            *(cti(iphbv)+ctk(iphbv))/2._10
                dvert(1)=dvert(1)/2._10
                dvert(2)=dvert(2)/2._10
              endif
              if (iphbv.le.2) then
                vert(1)=vert(1)+dvert(1)
                vert(2)=vert(2)+dvert(2)
              else
                vert(1)=vert(1)-dvert(1)
                vert(2)=vert(2)-dvert(2)
              endif
C             print *,'vertex',i,k,' part ',iphbv,
C    *          ' psi =',psi(iphbv)
C             print *,'     dvert =',dvert(1),dvert(2),
C    *          ' vert =',vert(1),vert(2)
C             print *,' cot(',psi(iphbv),') =',1._10/tan(psi(iphbv)),
C    *          ' should =',ctpsi(iphbv)
C             print *,' tan(',psi(iphbv),'/2) =',tan(psi(iphbv)/2._10),
C    *          ' should =',t(iphbv)
C             print *,' cot(th_i) =',cti(iphbv),
C    *          ' cot(th_k) =',ctk(iphbv)
            enddo
          enddo
c        do another segment
        goto 220
      endif
  280 continue
c--------finish off
c        check angle is between 0 and 2*pi
      p=bound(1)/si/TWOPI
c     print *,rp(1,i),rp(2,i),rp(3,i),'angle/(2*pi)=',p
      if (p.lt.0._10) then
        print *,'*** from gphbv: angle/(2*pi)=',p,
     *    '  should be .ge. 0'
        warn=.true.
      elseif (p.gt.1._10) then
        if (bound(1)/si.le.TWOPI+bndtol) then
          continue
        else
          print *,'*** from gphbv: angle/(2*pi)=',p,
     *      '  should be .le. 1'
          warn=.true.
        endif
      endif
      if (warn) then
        print *,'boundary',i
        write (*,'(a3,a20,4a24)')
     *    ' ','x','y','z',
c    *    'r',
     *    '1-c'
        do j=1,np
          write (*,'(i3,5g24.16)')
     *      j,(rp(i,j),i=1,3),
c    *      sqrt(rp(1,j)**2+rp(2,j)**2+rp(3,j)**2),
     *      cm(j)
        enddo
      endif
  410 continue
C     print *,'final bound =',bound(1),bound(2),
C    *  ' vert =',vert(1),vert(2)
C     print *,'....................'
      return
c
  420 continue
      print *,'*** from gphbv: total failure at tol =',tol
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
C     print *,'....................'
      return
c
      end
c
