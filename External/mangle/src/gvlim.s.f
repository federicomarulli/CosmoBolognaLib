c-----------------------------------------------------------------------
c © A J S Hamilton 2001
c-----------------------------------------------------------------------
      subroutine gvlim(vmin,vmax,cmvmin,cmvmax,cmpmin,cmpmax,
     *  ipv,gp,ev,nvmax,nv,nev,nev0,
     *  rp,cm,np,rpi,vcirc,tol,phi,iord,wk,iwk,ldegen)
      integer nvmax,ipv(nvmax),np,gp(np),ev(nvmax),nv,nev,nev0,vcirc
      logical ldegen
      real*10 vmin(3,nvmax),vmax(3,nvmax),cmvmin(nvmax),cmvmax(nvmax),
     *  cmpmin(np),cmpmax(np),rp(3,np),cm(np),rpi(3),tol
c        work arrays (could be made automatic if compiler supports it)
      integer iord(2*np),iwk(nvmax,4)
      real*10 phi(2,np),wk(nvmax)
c
c        parameters
      include 'pi.par'
      real*10 TWOPI
      parameter (TWOPI=2._10*PI)
c        intrinsics
      intrinsic abs
c        externals
      integer gsegij,gzeroar
      real*10 cmijf
c        data variables
      real*10 big
c     real*10 dphmin
c        local variables
      integer i,ii,ik,iphi,iseg,iv,j,jm,jml,jmu,jp,jpl,jpu,km,kp,ni,scmi
C     logical warn
      real*10 cmi,cmik,dph,ikchk,ph,phin,phif,phm,phmax,phmin,php,
     *  si,sik,tolin,xi(3),xv,yi(3),yv,zv
c *
c * Lifted mostly from gvert and gcmlim.
c *
c * Points nearest to and farthest from unit direction rpi
c * on each edge of polygon defined by
c *    1 - r.rp(i) < cm(i)  (if cm(i).ge.0)
c *    1 - r.rp(i) > -cm(i)  (if cm(i).lt.0)
c * for i=1,np where rp(i) are unit directions.
c *
c * The edges are ordered right-handedly about the polygon
c * in the same order as given by gvert.
c * The array ev records the end index of each connected
c * sequences of vertices.
c * Intersecting boundaries come 1st, then non-intersecting boundaries.
c * If the polygon is connected and simply-connected, the usual case,
c * then ev(1) = nv on output, and ev(i) = 0 for i >= 2.
c * If the polygon is not connected and simply-connected,
c * then the number of non-zero elements of ev will equal the number
c * of distinct connected boundaries of the polygon.
c *
c  Input: nvmax = guess at maximum nv;
c                 must be >= np;
c                 recommend:
c               = np        if 0 <= np <= 2
c               = np*(np-1) if 2 <= np <= 4
c               = 6*(np-2)  if np >= 4
c         rp(3,i),i=1,np
c         cm(i),i=1,np
c         np
c         rpi(3) = unit vector.
c         vcirc = 1 to return nearest and farthest points also for
c                 non-intersecting boundary circles.
c               = 0 otherwise.
c Output: vmin(3,i),vmax(3,i),i=1,nv are unit vectors which are the
c            nearest and farthest points on each edge from rpi.
c         cmvmin, cmvmax = 1 - cos th between nearest and farthest
c                points vmin, vmax on each edge from rpi.
c         cmpmin, cmpmax = 1 - cos th between nearest and farthest
c                points on each circle from rpi.
c         ipv(iv),iv=1,nv is the circle number of the iv'th edge,
c                  that is,points ve(*,ive,iv) lie on circle number ipv.
c         gp(i),i=1,np = which group circle i belongs to:
c                circle i and circle j are friends if they intersect,
c                and friends of friends are friends.
c         ev(i) = end index of each connected sequence of edges;
c                 the number of non-zero elements of ev is the number
c                 nev of connected boundaries of the polygon,
c                 and the last nev0 of these are non-intersecting.
c         nv = number of points;
c              if this exceeds nvmax, then you should call gvlim again
c              with a larger nvmax.
c         nev = number of connected sequences of vertices,
c               including those in non-intersecting circles.
c         nev0 = number of non-intersecting boundary circles if vcirc=1;
c              = 0 if vcirc=0;
c                non-intersecting boundaries come last in ve,
c                and the last nev0 entries of ev(i) refer
c                to non-intersecting boundaries.
c         ldegen = .true. means there's a problem with multiply
c                  intersecting boundary.
c Input/Output: tol
c Work arrays: phi and iord should be dimensioned at least 2*np.
c              iwk should be dimensioned at least 4*nvmax.
c              wk should be dimensioned at least nvmax.
c
c     data dphmin /1.e-8_10/
      data big /1.e6_10/
c
c        input tolerance to multiple intersections
      tolin=tol
c        come here with modified tolerance
  100 continue
c        initialise error flag to no error
      ldegen=.false.
C     warn=.false.
c        zero number of vertices
      nv=0
      nev=0
      nev0=0
c        initialise ev to zero
      do iv=1,nvmax
        ev(iv)=0
      enddo
c        initially each circle is its own group
      do i=1,np
        gp(i)=i
      enddo
c        check for zero area because one circle is null
      if (gzeroar(cm,np).eq.0) goto 410
c        error check on evaluation of vertex terms
      ikchk=0._10
c        initialise iwk to inadmissible value
      do iv=1,nvmax
        iwk(iv,2)=-1
      enddo
      do iv=1,nvmax
        iwk(iv,3)=-1
      enddo
c--------identify boundary segments around each circle i in turn
      do 280 i=1,np
c        cm(i).ge.2 means include whole sphere, which is no constraint
        if (cm(i).ge.2._10) then
          cmpmin(i)=2._10
          cmpmax(i)=0._10
          goto 280
        endif
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
        cmik=cmijf(rpi,rp(1,i))
c        sik = sin th(ik)
        sik=sqrt(cmik*(2._10-cmik))
c........minimum and maximum cm on circle
        cmpmin(i)=cmi+cmik-cmi*cmik-si*sik
        cmpmax(i)=cmi+cmik-cmi*cmik+si*sik
c........cartesian axes with z-axis along rp(i)
        call gaxisi(rp(1,i),xi,yi)
c........azimuthal angle closest to vector rpi
        xv=xi(1)*rpi(1)+xi(2)*rpi(2)+xi(3)*rpi(3)
        yv=yi(1)*rpi(1)+yi(2)*rpi(2)+yi(3)*rpi(3)
        phin=atan2(yv,xv)
        if (phin.ge.0._10) then
          phif=phin-PI
        else
          phif=phin+PI
        endif
c........angles phi about z-axis rp(i) of intersection of i & j circles
        call gphij(rp,cm,np,i,rp(1,i),scmi,cmi,xi,yi,big,tol,ni,phi)
c        i circle lies outside polygon
        if (ni.eq.-1) goto 280
c        area of polygon is zero
        if (ni.eq.-2) then
          nv=0
          nev=0
          nev0=0
          goto 410
        endif
c........i circle has no intersections
        if (ni.eq.0) then
c        want near and far points for non-intersecting boundary circles
          if (vcirc.eq.1) then
c        introduce pretend circle 0
            do j=1,2
              if (j.eq.1) then
                ii=0
                km=i
                kp=i
                phm=0._10
                php=PI
              elseif (j.eq.2) then
                ii=i
                km=0
                kp=0
                phm=PI
                php=TWOPI
              endif
              ph=(phm+php)/2._10
              dph=php-phm
              nv=nv+1
              if (nv.le.nvmax) then
c        edge index of this vertex
                ipv(nv)=i
c        phase phin to central point ph
                iphi=nint((phin-ph)/TWOPI)
                phin=phin-iphi*TWOPI
                if (phm.le.phin.and.phin.le.php) then
                  phmin=0._10
                elseif (phm.gt.phin) then
                  phmin=phm-phin
                elseif (phin.gt.php) then
                  phmin=php-phin
                endif
c        phase phif to central point ph
                iphi=nint((phif-ph)/TWOPI)
                phif=phif-iphi*TWOPI
                if (phm.le.phif.and.phif.le.php) then
                  phmax=0._10
                elseif (phm.gt.phif) then
                  phmax=phm-phif
                elseif (phif.gt.php) then
                  phmax=php-phif
                endif
c        nearest point on edge
                xv=si*cos(phmin+phin)
                yv=si*sin(phmin+phin)
                zv=1._10-cmi
                vmin(1,nv)=zv*rp(1,i)+xv*xi(1)+yv*yi(1)
                vmin(2,nv)=zv*rp(2,i)+xv*xi(2)+yv*yi(2)
                vmin(3,nv)=zv*rp(3,i)+xv*xi(3)+yv*yi(3)
c        farthest point on edge
                xv=si*cos(phmax+phif)
                yv=si*sin(phmax+phif)
                zv=1._10-cmi
                vmax(1,nv)=zv*rp(1,i)+xv*xi(1)+yv*yi(1)
                vmax(2,nv)=zv*rp(2,i)+xv*xi(2)+yv*yi(2)
                vmax(3,nv)=zv*rp(3,i)+xv*xi(3)+yv*yi(3)
c        minimum, maximum cm
                cmvmin(nv)=cmijf(rpi,vmin(1,nv))
                cmvmax(nv)=cmijf(rpi,vmax(1,nv))
              endif
c        record endpoints of edge
              call gvtrail(scmi,np,ii,km,kp,iwk,iwk(1,2),nvmax,nv,
     *          ik,ikchk)
            enddo
            nev0=nev0+1
          endif
c........i circle has intersections
        elseif (ni.gt.0) then
c........friends of i circle
          call ggpij(np,gp,i,big,phi)
c........find ordering of intersection angles around i circle
          call findbot(phi,2*np,iord,ni)
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
            if (iseg.eq.2) goto 280
c........segment is edge
c        warn about near multiple intersection
c           if (dph.lt.dphmin) then
c             print *,'*** warning from gvlim: near multiple intersectio
c    *n at',i,': edge',km,kp,' dph=',real(dph)
c             warn=.true.
c           endif
c           print *,'at',i,': edge',km,kp,
c    *        ' (',jm,' in',jml,jmu,',',jp,' in',jpl,jpu,' of',ni,')',
c    *        ' ph=',real(ph),' dph=',real(dph)
            nv=nv+1
            if (nv.le.nvmax) then
c        edge index of this vertex
              ipv(nv)=i
c        phase phin to central point ph
              iphi=nint((phin-ph)/TWOPI)
              phin=phin-iphi*TWOPI
              if (phm.le.phin.and.phin.le.php) then
                phmin=0._10
              elseif (phm.gt.phin) then
                phmin=phm-phin
              elseif (phin.gt.php) then
                phmin=php-phin
              endif
c        phase phif to central point ph
              iphi=nint((phif-ph)/TWOPI)
              phif=phif-iphi*TWOPI
              if (phm.le.phif.and.phif.le.php) then
                phmax=0._10
              elseif (phm.gt.phif) then
                phmax=phm-phif
              elseif (phif.gt.php) then
                phmax=php-phif
              endif
c        nearest point on edge
              xv=si*cos(phmin+phin)
              yv=si*sin(phmin+phin)
              zv=1._10-cmi
              vmin(1,nv)=zv*rp(1,i)+xv*xi(1)+yv*yi(1)
              vmin(2,nv)=zv*rp(2,i)+xv*xi(2)+yv*yi(2)
              vmin(3,nv)=zv*rp(3,i)+xv*xi(3)+yv*yi(3)
c        farthest point on edge
              xv=si*cos(phmax+phif)
              yv=si*sin(phmax+phif)
              zv=1._10-cmi
              vmax(1,nv)=zv*rp(1,i)+xv*xi(1)+yv*yi(1)
              vmax(2,nv)=zv*rp(2,i)+xv*xi(2)+yv*yi(2)
              vmax(3,nv)=zv*rp(3,i)+xv*xi(3)+yv*yi(3)
c        minimum, maximum cm
              cmvmin(nv)=cmijf(rpi,vmin(1,nv))
              cmvmax(nv)=cmijf(rpi,vmax(1,nv))
            endif
c        record endpoints of edge
            call gvtrail(scmi,np,i,km,kp,iwk,iwk(1,2),nvmax,nv,ik,ikchk)
c        do another segment
          goto 220
        endif
  280 continue
c--------check on whether ik endpoints matched ki endpoints
      if (ikchk.ne.0._10) then
c       print *,'*** from gvlim: ikchk=',ikchk,' should be 0'
C       warn=.true.
        call gtol(tol,tolin)
        goto 100
      endif
c--------order vertices right-handedly about polygon
      if (nv.gt.0.and.nv.le.nvmax) then
c........which group of friends each circle belongs to
        call ggp(np,gp)
c........find order of vertices around polygon
        call gvord(np,nv,ipv,gp,ev,nev,iwk,iwk(1,2),iwk(1,3),iwk(1,4))
c........reorder vertices
        call vpermi(ipv,nv,iwk(1,2),iwk(1,3))
        call vpermdd(vmin,3,nv,iwk(1,2),wk)
        call vpermdd(vmax,3,nv,iwk(1,2),wk)
        call vpermd(cmvmin,nv,iwk(1,2),wk)
        call vpermd(cmvmax,nv,iwk(1,2),wk)
      endif
c--------done
 410  continue
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
c     if (nv.gt.nvmax) then
c       print *,'*** from gvlim: number of vertices =',nv,
c    *    ' exceeds maximum',nvmax
c     endif
      return
c
  420 print *,'*** from gvlim: total failure at tol =',tol
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
      ldegen=.true.
      return
c
      end
c
