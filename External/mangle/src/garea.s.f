c-----------------------------------------------------------------------
c © A J S Hamilton 2001
c-----------------------------------------------------------------------
      subroutine garea(area,rp,cm,np,tol,verb,phi,iord,ldegen)
      integer np,verb
      logical ldegen
      real*10 area,rp(3,np),cm(np),tol
c        work arrays (could be automatic if compiler supports it)
      integer iord(2*np)
      real*10 phi(2,np)
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
c        local variables
      integer i,iarea,ik,iseg,j,jm,jml,jmu,jp,jpl,jpu,k,km,kp,l,
     *  nbd,nbd0m,nbd0p,ni,nmult,retry,scmi
C     logical warn
      logical whole
      real*10 bik,cmi,cmik,cmk,d,darea,dph,ikchk,ikran,
     *  ph,phm,php,psi,si,tolin,xi(3),yi(3)
c *
c * Area of surface of sphere of unit radius bounded by
c *    1 - r.rp(i) < cm(i)  (if cm(i).ge.0)
c *    1 - r.rp(i) > -cm(i)  (if cm(i).lt.0)
c * for i=1,np where rp(i) are unit directions.
c * See AJSH notes Multfn C115.
c * Cautions:
c * (1) This subroutine underestimates the area by 2*pi
c *     if 2*pi <= area < 4*pi
c *     and the area is bounded by more than one arc.
c * (2) This subroutine will usually work correctly when there are near
c *     multiple (.ge. 3) intersections of boundaries, but in rare
c *     instances it may fail.  If so, it should flag the failure with
c *     ldegen=.true.  This error condition should NOT be ignored.
c *
c  Input: rp(3,i),i=1,np
c         cm(i),i=1,np
c         np
c         verb
c Output: area
c         ldegen = .true. means there's a problem with multiply
c                  intersecting boundary.
c Input/Output: tol
c Work arrays: phi and iord should be dimensioned at least 2*np
c
c        set azimuthal angle of non-intersection to big
      data big /1.e6_10/
c        possible multiple intersection when dph < dphmin
      data dphmin /1.e-8_10/
c
c        input tolerance to multiple intersections
      tolin=tol
C     print *,'--------------------'
c        come here with modified tolerance to multiple intersections
  100 continue
C     write (*,'(a3,a20,4a24)')
C    *  ' ','x','y','z',
c    *  'r',
C    *  '1-c'
C     do j=1,np
C       write (*,'(i3,5g24.16)')
C    *    j,(rp(i,j),i=1,3),
c    *    sqrt(rp(1,j)**2+rp(2,j)**2+rp(3,j)**2),
C    *    cm(j)
C     enddo
c        initialise error flag to no error
      ldegen=.false.
C     warn=.false.
c        initialize count of near multiple intersections to zero
      nmult=0
c        zero area
      area=0._10
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
c        si = sin th(i)
        si=sqrt(cmi*(2._10-cmi))
c........construct cartesian axes with z-axis along rp(i)
        call gaxisi(rp(1,i),xi,yi)
c........angles phi about z-axis rp(i) of intersection of i & j circles
        call gphij(rp,cm,np,i,rp(1,i),scmi,cmi,xi,yi,big,tol,ni,phi)
c        i circle lies outside polygon
        if (ni.eq.-1) goto 280
c        area of polygon is zero
        if (ni.eq.-2) then
c        area can be non-zero from psi at multiple intersections
          area=0._10
          goto 410
        endif
c........i circle has no intersections
        if (ni.eq.0) then
          if (scmi.ge.0) then
            nbd0p=nbd0p+1
          else
            nbd0m=nbd0m+1
          endif
          darea=cm(i)*TWOPI
          area=area+darea
C         print *,'at',i,': full circle area +=',darea,' =',area
c........i circle has intersections
        elseif (ni.gt.0) then
c        find ordering of intersection angles around i circle
          call findbot(phi,2*np,iord,ni)
C         write (*,'("phi")')
C         write (*,'(i4,2g24.16)') (j,(phi(k,j),k=1,2),j=1,np)
c........contribution to area from each segment of i circle
          jpl=0
c        come here to do another segment
  220     continue
c........is segment edge of polygon?
            iseg=gsegij(rp,cm,np,0,0,i,rp(1,i),scmi,cmi,tol,ni,
     *        phi,iord,jml,jmu,jpl,jpu,1,jm,jp,km,kp,phm,php,ph,dph)
C           print *,'at',i,': iseg =',iseg
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
C             warn=.true.
C             print *,
C    *          '*** warning from garea: near multiple intersection at',
C    *          i,': edge',km,kp,' dph =',dph
            endif
c........segment satisfies conditions
            nbd=nbd+1
c        there's a contribution to area from the segment...
            if (scmi.lt.0) dph=-dph
            darea=cmi*dph-dph
            area=area+darea
C           print *,'at',i,': edge',km,kp,
C    *       ' (',jm,' in',jml,jmu,',',jp,' in',jpl,jpu,' of',ni,')'
C           print *,'dph =',dph,' area +=',darea,' =',area
c        ...and from the end points of the segment
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
C               print *,'    (intersect',i,k,')'
                goto 240
              endif
c        ikchk = ikchk + ikran, added as unsigned long long's
              call ikrandp(ikchk,ikran)
              cmk=abs(cm(k))
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
              if (phi(1,k).eq.phi(2,k)) then
                d=0._10
              else
                d=-(cmi-cmk)**2+cmik*(2._10*((cmi+cmk)-cmi*cmk)-cmik)
c        assert that circles at least touch
                if (d.lt.0._10) d=0._10
                d=sqrt(d)
              endif
              psi=atan2(d,bik)
              area=area-psi
C             print *,'     intersect',i,k,' area +=',-psi,' =',area
  240       continue
c        do another segment
          goto 220
        endif
  280 continue
c--------check on whether ik endpoints matched ki endpoints
      if (ikchk.ne.0._10) then
C       warn=.true.
C       print *,'*** from garea: at tol =',tol,
C    *    ', ikchk=',ikchk,' should be 0'
c       write (*,'(a3,a20,3a24)')
c    *    ' ','x','y','z','1-c'
c       write (*,'(i3,4g24.16)')
c    *    (j,(rp(i,j),i=1,3),
c    *    cm(j),j=1,np)
c        retry with modified tolerance
        call gtol(tol,tolin)
        goto 100
      elseif (tol.gt.0._10) then
C       print *,'... from garea: success at tol =',tol
      endif
c--------add/subtract 2*pi's to area
      retry=garpi(area,iarea,rp,cm,np,whole,nbd0m,nbd0p,nbd,nmult)
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
C     print *,'final area =',area
C     print *,'....................'
      return
c
  420 continue
c     if (verb.ge.1) then
        print *,'*** from garea: total failure at tol =',tol
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
c     endif
C     print *,'....................'
      ldegen=.true.
      return
c
      end
c
