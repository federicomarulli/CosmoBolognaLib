c-----------------------------------------------------------------------
c © A J S Hamilton 2001
c-----------------------------------------------------------------------
      subroutine gspher(area,bound,vert,w,lmax1,im,nw,rp,cm,np,npc,ibv,
     *  iphi,tol,phw,iord,v,ldegen)
      integer lmax1,im,nw,np,npc,ibv,iphi
      logical ldegen
      real*10 area,bound(2),vert(2),w(im,nw),rp(3,np),cm(np),tol
c        work arrays (could be automatic if compiler supports it)
      integer iord(2*np)
      real*10 phw(2,np),v(lmax1)
c
c        parameters
      include 'pi.par'
      real*10 TWOPI
      parameter (TWOPI=2._10*PI)
c        intrinsics
      intrinsic abs
c        externals
      integer garpi,gsegij,gzeroar
c        data variables
      real*10 big
      real*10 dphmin
c        local (automatic) variables
      integer i,iarea,ik,iseg,j,jm,jml,jmu,jp,jpl,jpu,k,km,kp,l,
     *  nbd,nbd0m,nbd0p,ni,nmult,retry,scmi
C     logical warn
      logical whole
      real*10 bik,ci,cmi,cmik,cmk,cti,ctk,ctpsi,
     *  d,darea,dbound(2),dph,dvert(2),ikchk,ikran,
     *  ph,phi,phii,phm,php,psi,psip,ri,rii,si,sk,sqrt4pi,t,tolin,
     *  xi(3),yi(3)
c *
c * Spherical transform of region W of sphere of unit radius bounded by
c *    1 - r.rp(i) < cm(i)  (if cm(i).ge.0)
c *    1 - r.rp(i) > -cm(i)  (if cm(i).lt.0)
c * for i=1,np where rp(i) are unit directions.
c * See AJSH notes Multfn C115.
c *
c * bound (boundary) and vert (vertex) are related to the coefficients
c * of the 1st & 2nd order terms of the power series expansion
c * in sin(th/2) of the correlation <WW> at angular separation th:
c *    <WW> = sum(lm) |W_lm|^2 P_l(cos th)
c *         = 2*pi*area - 4*bound*sin(th/2) + 2*vert*sin^2(th/2) + ...
c * In fact bound is just the total length of the boundary in radians,
c * and vert is the sum over vertex terms of 1-psi/tan(psi), where
c * psi is the exterior angle (=pi-interior angle) at the vertex.
c * The range of validity of the power series expansion in sin(th/2)
c * depends on the tortuousness of the boundary: the more convoluted the
c * boundary, and the sharper the vertex angles, the shorter the range.
c *
c * Oct 1992: next order term incorporated
c *    <WW> = 2*pi*area - 4*bound(1)*sin(th/2) + 2*vert(1)*sin^2(th/2)
c *           + [2/3*bound(2) + 8/9*vert(2)]*sin^3(th/2) + ...
c *
c * The option ibv modifies the evaluation of bound and vert,
c * and determines the sign of the returned harmonics.
c * ibv = 0: standard option, for <WW>.
c *          The code returns the harmonics of W in the array w.
c * ibv = 1: to take the cross-correlation between a region W1 and
c *          its intersection W12 with region W2.
c *          Here <W12 W1> = <W12^2> + <W12(W1-W12)>
c *          and this is what the code evaluates bound and vert for.
c *          The code returns the harmonics of W12 in the array w.
c * ibv = 2: to take region W1 less its intersection W12 with region W2.
c *          Here <WW> = <(W1-W12)(W1-W12)>
c *                    = <W1^2> - <W12^2> - 2<W12(W1-W12)>
c *          The code evaluates bound and vert for
c *                               <W12^2> + 2<W12(W1-W12)>.
c *          The code returns MINUS the harmonics of W12 in the array w.
c * ibv = 3: to take the union of regions W1 and W2.
c *          Here <WW> = <(W1+W2-W12)(W1+W2-W12)>
c *                    = <W1^2> + <W2^2> - <W12^2> + 2<(W1-W12)(W2-W12)>
c *          The code evaluates bound and vert for
c *                               <W12^2> - 2<(W1-W12)(W2-W12)>.
c *          The code returns MINUS the harmonics of W12 in the array w.
c * The correct calling procedure in each case is exampled below.
c *
c * Complement of a region:
c *    area -> 4*pi-area,  w(1,1) -> sqrt(4*pi)-w(1,1)
c *    bound -> bound
c *    vert -> vert
c *    w -> -w   except monopole term w(1,1) as above
c *
c * To get intersection W12 of regions W1 & W2,
c * but evaluate bound & vert for cross-correlation between W12 and W1:
c *    put the np2 constraints of region W2 in 1 to npc
c *      & the np1 constraints of region W1 in npc+1 to np
c *    call gspher(area,bound,vert,w,parameters of W1,npc=np2,1,...)
c *
c * To get region W1 less its intersection with region W2:
c *    call gspher(area,bound,vert,w,parameters of W1,npc=0,0,...)
c *    put the np2 constraints of region W2 in 1 to npc
c *      & the np1 constraints of region W1 in npc+1 to np
c *    call gspher(darea,dbound,dvert,w,params of W1 & W2,npc=np2,2,...)
c *    area=area-darea
c *    bound=bound-dbound
c *    vert=vert-dvert
c *
c * To get the union of two regions W1 & W2:
c *    call gspher(area,bound,vert,w,parameters of W1,npc=0,0,...)
c *    call gspher(darea,dbound,dvert,w,parameters of W2,npc=0,0,...)
c *    area=area+darea
c *    bound=bound+dbound
c *    vert=dvert+dvert
c *    put the np2 constraints of region W2 in 1 to npc
c *      & the np1 constraints of region W1 in npc+1 to np
c *    call gspher(darea,dbound,dvert,w,params of W1 & W2,npc=np2,3,...)
c *    area=area-darea
c *    bound=bound-dbound
c *    vert=dvert-dvert
c *
c * Cautions:
c * (1) This subroutine underestimates the area (monopole harmonic)
c *     by 2*pi (sqrt(pi)) if 2*pi <= area < 4*pi
c *     and the area is bounded by more than one arc.
c * (2) This subroutine will usually work correctly when there are near
c *     multiple (.ge. 3) intersections of boundaries, but in rare
c *     instances it may fail.  If so, it should flag the failure with
c *     ldegen=.true.  This error condition should NOT be ignored.
c * (3) There is a mathematical ambiguity in the boundary term bound
c *     whenever two (or more) boundaries coincide.  To resolve the
c *     ambiguity correctly, perturb the boundaries slightly.
c *     Likewise there is a mathematical ambiguity in the vertex term
c *     vert whenever there are multiple (.ge. 3) intersections.
c *     Again, to resolve the ambiguity correctly, perturb the
c *     boundaries slightly.
c *     It should be noted that if these problems exist, then the series
c *     expansion of <WW>, whose coefficients involve bound and vert,
c *     probably breaks down already at tiny values of the separation
c *     angle th.
c *
c  Input: lmax1 = lmax+1 where lmax is maximum desired l of transform.
c         im = 1 means compute only real part of harmonics;
c              2 means compute both real and imaginary parts.
c              Note harmonics are real if region possesses reflection
c              symmetry through plane defined by z-axis and direction
c              rp(iphi).
c         nw = [(lmax+1)*(lmax+2)]/2
c         rp(3,i),i=1,np = x, y, z coordinates of a set of unit
c              directions defining the region.
c              It is assumed without checking that
c              rp(1)^2 + rp(2)^2 + rp(3)^2 = 1 .
c         cm(i),i=1,np = set of 1-cos's defining the region.
c         np = number of directions.
c         npc = part number of directions, used in conjunction with ibv;
c               see above for more details;
c               npc is irrelevant if ibv = 0.
c         ibv = 0 to 3 controls evaluation of bound & vert;
c               see above for more details.
c         iphi > 0 means compute harmonics in frame of reference where
c                  y-axis is along z x rp(iphi);
c              = 0 means use the input frame of reference.
c Output: area = area of region in steradians.
c         bound = length of boundary of region in radians if ibv=0,
c                 or as explained above if ibv>0.
c         vert = sum over vertices of 1-psi/tan(psi) if ibv=0,
c                where psi is exterior angle (=pi-interior angle)
c                at vertex, or as explained above if ibv>0.
c         ldegen = .true. signals an error: the code dealt incorrectly
c                  with a multiply intersecting boundary.
c Input/Output: w(i,lm) = spherical transform, dimensioned w(im,nw)
c            w(i,lm), i=1,im, lm=l*(l+1)/2+m+1, l=0,lmax, m=0,l;
c            w(1,lm) is real part, w(2,lm) is imaginary part (if im=2).
c            Note w(l,-m)=(-)**m*[Complex conjugate of w(l,m)], just as
c                 Y(l,-m)=(-)**m*[Complex conjugate of Y(l,m)].
c            tol
c Work arrays: phw and iord should be dimensioned at least 2*np.
c              v should be dimensioned at least lmax1.
c
c        set azimuthal angle of non-intersection to big
      data big /1.e6_10/
c        possible multiple intersection when dph < dphmin
      data dphmin /1.e-8_10/
c
c        input tolerance to multiple intersections
      tolin=tol
C     print *,'--------------------'
c        come here with modified tolerance
  100 continue
c        initialise error flag to no error
      ldegen=.false.
C     warn=.false.
c        zero stuff
      area=0._10
      bound(1)=0._10
      bound(2)=0._10
      vert(1)=0._10
      vert(2)=0._10
      do j=1,nw
        do i=1,im
          w(i,j)=0._10
        enddo
      enddo
c        check for zero area because one circle is null
      if (gzeroar(cm,np).eq.0) goto 410
c        no constraints at all will mean area is whole sphere
      whole=.true.
c        number of intersecting arc segments bounding area
      nbd=0
c        number of non-intersecting circles bounding area
      nbd0m=0
      nbd0p=0
c        error check on evaluation of vertex terms
      ikchk=0._10
c        area=sqrt(4pi)*monopole
      sqrt4pi=sqrt(4._10*PI)
c        harmonics defined so point iphi is at zero azimuthal angle
      rii=0._10
      if (iphi.ge.1) rii=sqrt(rp(1,iphi)**2+rp(2,iphi)**2)
      phii=0._10
      if (rii.gt.0._10) phii=atan2(rp(2,iphi),rp(1,iphi))
c--------identify boundary segments around each circle i in turn
      do 280 i=1,np
c        cm(i).ge.2 means include whole sphere, which is no constraint
        if (cm(i).ge.2._10) goto 280
c        there is a constraint, so area is not whole sphere
        whole=.false.
c        scmi * cmi = 1-cos th(i)
        if (cm(i).ge.0._10) then
          scmi=1
        else
          scmi=-1
        endif
        cmi=abs(cm(i))
c        ci = cos th(i)
        ci=1._10-cmi
c        si = sin th(i)
        si=sqrt(cmi*(2._10-cmi))
c........ri, phi, rp(3,i) are cylindrical coordinates of rp(i)
        ri=sqrt(rp(1,i)**2+rp(2,i)**2)
        if (ri.eq.0._10.or.i.eq.iphi) then
          phi=0._10
        else
          phi=atan2(rp(2,i),rp(1,i))-phii
        endif
c........construct cartesian axes with z-axis along rp(i)
c The direction of yi is important here,
c unlike some other subroutines (gphi, garea, gphbv, gvlim, gvphi)
c where yi can point in any abitrary direction.
c        set yi in direction z x rp(i)
        if (ri.gt.0._10) then
          yi(1)=-rp(2,i)/ri
          yi(2)=rp(1,i)/ri
          yi(3)=0._10
c        if rp(i) is along z-axis, set yi in direction z x rp(iphi)
        elseif (rii.gt.0._10) then
          yi(1)=-rp(2,iphi)/rii
          yi(2)=rp(1,iphi)/rii
          yi(3)=0._10
c        if rp(iphi) is also along z-axis, set yi along y-axis
        elseif (ri.eq.0._10.and.rii.eq.0._10) then
          yi(1)=0._10
          yi(2)=1._10
          yi(3)=0._10
        endif
c        xi in direction yi x rp(i)
        xi(1)=yi(2)*rp(3,i)-yi(3)*rp(2,i)
        xi(2)=yi(3)*rp(1,i)-yi(1)*rp(3,i)
        xi(3)=yi(1)*rp(2,i)-yi(2)*rp(1,i)
c........angles phi about z-axis rp(i) of intersection of i & j circles
        call gphij(rp,cm,np,i,rp(1,i),scmi,cmi,xi,yi,big,tol,ni,phw)
c        i circle lies outside polygon
        if (ni.eq.-1) goto 280
c        area of polygon is zero
        if (ni.eq.-2) then
c        area can be non-zero from psi at multiple intersections
          area=0._10
          bound(1)=0._10
          bound(2)=0._10
          vert(1)=0._10
          vert(2)=0._10
          do l=1,nw
            do k=1,im
              w(k,l)=0._10
            enddo
          enddo
          goto 410
        endif
c........i circle has no intersections
        if (ni.eq.0) then
          if (scmi.ge.0) then
            nbd0p=nbd0p+1
          else
            nbd0m=nbd0m+1
          endif
          dph=TWOPI
          if (scmi.lt.0) dph=-dph
          ph=0._10
c        increment area
          darea=cmi*dph
          area=area+darea
c        bound(1) term is length of boundary
          dbound(1)=si*abs(dph)
          dbound(2)=(1._10/si-2._10*si)*abs(dph)
c        standard
          if (ibv.eq.0
c        cross
     *      .or.(ibv.eq.1.and.i.gt.npc)
c        intersection
     *      .or.(ibv.eq.2.and.i.gt.npc)
c        union
     *      .or.(ibv.eq.3)) then
            bound(1)=bound(1)+dbound(1)
            bound(2)=bound(2)+dbound(2)
c        intersection
          elseif (ibv.eq.2.and.i.le.npc) then
            bound(1)=bound(1)-dbound(1)
            bound(2)=bound(2)-dbound(2)
          endif
C         print *,'at',i,': full circle area +=',darea,' =',area
C         print *,'dbound =',dbound(1),dbound(2),
C    *      ' bound =',bound(1),bound(2)
c        increment spherical transform
          if (ibv.ge.2) dph=-dph
          call wlm(w,lmax1,im,nw,ri,phi,0,rp(3,i),ci,si,ph,dph,v)
c........i circle has intersections
        elseif (ni.gt.0) then
c        find ordering of intersection angles around i circle
          call findbot(phw,2*np,iord,ni)
c........contribution from each segment of i circle
          jpl=0
c        come here to do another segment
  220     continue
c........is segment edge of polygon?
            iseg=gsegij(rp,cm,np,0,0,i,rp(1,i),scmi,cmi,tol,ni,
     *        phw,iord,jml,jmu,jpl,jpu,1,jm,jp,km,kp,phm,php,ph,dph)
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
c             warn=.true.
c             print *,
c    *          '*** warning from gspher: near multiple intersection at'
c    *          ,i,': segment',km,kp,' dph=',dph
            endif
c........segment satisfies conditions
            nbd=nbd+1
c        contribution to area from the boundary segment
            if (scmi.lt.0) dph=-dph
            darea=cmi*dph-dph
            area=area+darea
c        bound(1) term is length of boundary
            dbound(1)=si*abs(dph)
            dbound(2)=(1._10/si-2._10*si)*abs(dph)
c        standard
            if (ibv.eq.0
c        cross
     *        .or.(ibv.eq.1.and.i.gt.npc)
c        intersection
     *        .or.(ibv.eq.2.and.i.gt.npc)
c        union
     *        .or.(ibv.eq.3)) then
              bound(1)=bound(1)+dbound(1)
              bound(2)=bound(2)+dbound(2)
c        intersection
            elseif (ibv.eq.2.and.i.le.npc) then
              bound(1)=bound(1)-dbound(1)
              bound(2)=bound(2)-dbound(2)
            endif
C           print *,'at',i,': edge',km,kp,
C    *       ' (',jm,' in',jml,jmu,',',jp,' in',jpl,jpu,' of',ni,')'
C           print *,'dph =',dph,' area +=',darea,' =',area
C           print *,'dbound =',dbound(1),dbound(2),
C    *        ' bound =',bound(1),bound(2)
c        contribution to area and vert from end points of the segment
            do 240 l=1,2
c        end point is intersection of i circle with k circle
              if (l.eq.1) then
                k=km
              elseif (l.eq.2) then
                k=kp
              endif
c        only do ik intersection once,
c        but check both ik and ki intersections made it here
c        from segment k to segment i right-handedly through vertex
              if ((l.eq.1.and.scmi.ge.0).or.(l.eq.2.and.scmi.lt.0)) then
                ik=k+np*i
                if (k.lt.i) ik=ik+1
c        from segment i to segment k right-handedly through vertex
              else
                ik=i+np*k
                if (i.lt.k) ik=ik+1
              endif
c        pseudo-random number from ik
              call ikrand(ik,ikran)
              if (i.gt.k) then
c        ikchk = ikchk - ikran, subtracted as unsigned long long's
                call ikrandm(ikchk,ikran)
                goto 240
              endif
c        ikchk = ikchk + ikran, added as unsigned long long's
              call ikrandp(ikchk,ikran)
              cmk=abs(cm(k))
              sk=sqrt(cmk*(2._10-cmk))
c        cmik = 1-cos th(ik)
              cmik=((rp(1,i)-rp(1,k))**2+(rp(2,i)-rp(2,k))**2
     *          +(rp(3,i)-rp(3,k))**2)/2._10
c        bik = cik-ci*ck
c        d = 1-ci^2-ck^2-cik^2+2*ci*ck*cik
c        cos psi = bik/(si*sk)
c        sin psi = sqrt(d)/(si*sk)
c        psi = atan(sqrt(d)/bik) is exterior angle at intersection
              bik=(cmi+cmk)-cmi*cmk-cmik
              if ((scmi.ge.0.and.cm(k).lt.0._10)
     *          .or.(scmi.lt.0.and.cm(k).ge.0._10)) bik=-bik
c        i and k circles kiss
              if (phw(1,k).eq.phw(2,k)) then
                d=0._10
              else
                d=-(cmi-cmk)**2+cmik*(2._10*((cmi+cmk)-cmi*cmk)-cmik)
c        assert that circles at least touch
                if (d.lt.0._10) d=0._10
                d=sqrt(d)
              endif
              ctpsi=bik/d
              psi=atan2(d,bik)
c        increment area
              area=area-psi
c        t=tan psi/2
              if (bik.gt.0._10) then
                t=d/(bik+sqrt(bik**2+d**2))
              elseif (bik.lt.0._10) then
                t=(-bik+sqrt(bik**2+d**2))/d
              elseif (bik.eq.0._10) then
                t=1._10
              endif
c        cti = cot th(i)
              cti=(1._10-cmi)/si
              if (scmi.lt.0) cti=-cti
c        ctk = cot th(k)
              ctk=(1._10-cmk)/sk
              if (cm(k).lt.0._10) ctk=-ctk
c        standard
              if (ibv.eq.0
     *          .or.(ibv.eq.1.and.i.gt.npc.and.k.gt.npc)
     *          .or.(ibv.eq.2.and.i.gt.npc.and.k.gt.npc)
     *          .or.(ibv.eq.3.and.((i.gt.npc.and.k.gt.npc)
     *                         .or.(i.le.npc.and.k.le.npc)))) then
                if (psi.eq.0._10) then
                  dvert(1)=0._10
                  dvert(2)=0._10
                else
                  dvert(1)=1._10-psi*ctpsi
                  dvert(2)=t*(3._10+t**2)*(cti+ctk)/2._10
                endif
                vert(1)=vert(1)+dvert(1)
                vert(2)=vert(2)+dvert(2)
c        cross
              elseif (ibv.eq.1) then
                if (i.le.npc.and.k.le.npc) then
                  continue
                else
                  dvert(1)=PI/2._10*ctpsi
                  dvert(2)=t*(3._10+t**2)*(cti+ctk)/2._10
                  vert(1)=vert(1)-dvert(1)
                  vert(2)=vert(2)+dvert(2)/2._10
                  t=1._10/t
                  dvert(2)=t*(3._10+t**2)*(cti-ctk)/2._10
                  if (i.gt.npc) then
                    vert(2)=vert(2)-dvert(2)/2._10
                  elseif (k.gt.npc) then
                    vert(2)=vert(2)+dvert(2)/2._10
                  endif
                endif
c        intersection
              elseif (ibv.eq.2) then
                if (i.le.npc.and.k.le.npc) then
                  if (psi.eq.0._10) then
                    dvert(1)=0._10
                    dvert(2)=0._10
                  else
                    dvert(1)=1._10-psi*ctpsi
                    dvert(2)=t*(3._10+t**2)*(cti+ctk)/2._10
                  endif
                  vert(1)=vert(1)-dvert(1)
                  vert(2)=vert(2)-dvert(2)
                else
c        psip = pi - psi
                  psip=atan2(d,-bik)
                  if (psip.eq.0._10) then
                    dvert(1)=0._10
                    dvert(2)=0._10
                  else
                    dvert(1)=1._10+psip*ctpsi
                    t=1._10/t
                    dvert(2)=t*(3._10+t**2)*(cti-ctk)/2._10
                  endif
                  vert(1)=vert(1)-dvert(1)
                  if (i.gt.npc) then
                    vert(2)=vert(2)-dvert(2)
                  elseif (k.gt.npc) then
                    vert(2)=vert(2)+dvert(2)
                  endif
                endif
c        union
              elseif (ibv.eq.3) then
                if (psi.eq.0._10) then
                  dvert(1)=0._10
                  dvert(2)=0._10
                else
                  dvert(1)=1._10-psi*ctpsi
                  dvert(2)=t*(3._10+t**2)*(cti+ctk)/2._10
                endif
                vert(1)=vert(1)-dvert(1)
                vert(2)=vert(2)+dvert(2)
              endif
C             print *,'     intersect',i,k,' area +=',-psi,' =',area
C             print *,'     dvert =',dvert(1),dvert(2),
C    *          ' vert =',vert(1),vert(2)
C             print *,' cot(',psi,') =',1._10/tan(psi),
C    *          ' should =',ctpsi
C             print *,' tan(',psi,'/2) =',tan(psi/2._10),
C    *          ' should =',t 
C             print *,' cot(th_i) =',cti,' cot(th_k) =',ctk
c        peculiar monopole term
              if (ibv.ge.2) psi=-psi
              w(1,1)=w(1,1)-psi/sqrt4pi
  240       continue
c        increment spherical transform
            if (ibv.ge.2) dph=-dph
            call wlm(w,lmax1,im,nw,ri,phi,0,rp(3,i),ci,si,ph,dph,v)
c        do another segment
          goto 220
        endif
  280 continue
c--------check on whether ik endpoints matched ki endpoints
      if (ikchk.ne.0._10) then
C       warn=.true.
c       print *,'*** from gspher: at tol =',tol,
c    *    ', ikchk=',ikchk,' should be 0'
c        retry with enlarged tolerance
        call gtol(tol,tolin)
        goto 100
      endif
c--------add/subtract 2*pi's to area
      retry=garpi(area,iarea,rp,cm,np,whole,nbd0m,nbd0p,nbd,nmult)
c        adjust monopole harmonic by corresponding 2*pi/sqrt(4*pi)
      if (ibv.eq.0.or.ibv.eq.1) then
        i=nint((w(1,1)*sqrt4pi-area)/TWOPI)
      elseif (ibv.eq.2.or.ibv.eq.3) then
        i=nint((w(1,1)*sqrt4pi+area)/TWOPI)
      endif
      w(1,1)=w(1,1)-i*TWOPI/sqrt4pi
c        retry with modified tolerance
      if (retry.eq.1) then
C       warn=.true.
        call gtol(tol,tolin)
        goto 100
      endif
c--------done
  410 continue
C     if (warn) then
C       write (*,'(a3,a20,4a24)')
C    *    ' ','x','y','z',
c    *    'r',
C    *    '1-c'
C       do j=1,np
C         write (*,'(i3,5g24.16)')
C    *      j,(rp(i,j),i=1,3),
c    *      sqrt(rp(1,j)**2+rp(2,j)**2+rp(3,j)**2),
C    *      cm(j)
C       enddo
C     endif
C     print *,'final area =',area,' bound =',bound(1),bound(2),
C    *  ' vert =',vert(1),vert(2)
C     print *,'....................'
      return
c
  420 print *,'*** from gspher: total failure at tol =',tol
      ldegen=.true.
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
