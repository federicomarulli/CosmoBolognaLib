c-----------------------------------------------------------------------
c © A J S Hamilton 2001
c-----------------------------------------------------------------------
      subroutine gvert(ve,angle,ipv,gp,ev,nvmax,nv,per,nve,nev,nev0,
     *  rp,cm,np,vcirc,tol,phi,iord,wk,iwk,ldegen)
      integer nvmax,ipv(nvmax),np,gp(np),ev(nvmax),nv,per,nve,nev,nev0,
     *  vcirc
      logical ldegen
      real*10 ve(3,nve,nvmax),angle(nvmax),rp(3,np),cm(np),tol
c        work arrays (could be automatic if compiler supports it)
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
c        data variables
      real*10 big
c     real*10 dphmin
c        local (automatic) variables
      integer i,ii,ik,iseg,iv,ive,j,jm,jml,jmu,jp,jpl,jpu,km,kp,
     *  mve,ni,scmi
C     logical warn
      real*10 amve,cmi,dph,ikchk,ph,phm,php
      real*10 si,tolin,xi(3),xv,yi(3),yv,zv
c *
c * Vertices, plus points on edges, of the polygon defined by
c *    1 - r.rp(i) < cm(i)  (if cm(i).ge.0)
c *    1 - r.rp(i) > -cm(i)  (if cm(i).lt.0)
c * for i=1,np where rp(i) are unit directions.
c *
c * Since each pair of circles can intersect at at most 2 points,
c * an upper limit on the number of possible vertices is np*(np-1),
c * twice the number of distinct pairs of np objects.
c * I can find configurations which exhaust this bound for np <= 4.
c * For np >= 4, I can find configurations with 6*(np-2) vertices,
c * namely, form two triangles from the intersection of 3 caps
c * (producing 6 vertices),
c * then for each additional cap, inscribe it as an excluding cap
c * inside one of the triangles, so as to transform that one
c * triangle (with 3 vertices) into 3 triangles (with 9 vertices).
c *
c * gp(i) is the group of circles to which circle i belongs.
c * Two circles are friends, belonging to the same group,
c * if they intersect, and friends of friends are friends.
c * Theorem: A polygon is not simply-connected
c * if and only if its boundaries belong to distinct groups.
c * The group number is the smallest circle number of circles
c * within the group.
c *
c * Connected sequences of vertices are ordered by group,
c * with the lowest group numbers coming first, except that
c * intersecting boundaries come 1st, then non-intersecting boundaries.
c * Each connected sequence of vertices is ordered right-handedly
c * about the polygon, with the vertices belonging to the lowest
c * circle numbers coming first.
c * The array ev records the end index of each connected
c * sequences of vertices.
c * If the polygon is connected and simply-connected, the usual case,
c * then ev(1) = nv on output, and ev(i) = 0 for i >= 2.
c * If the polygon is not connected and simply-connected,
c * then the number of non-zero elements of ev will equal the
c * number of distinct connected boundaries of the polygon.
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
c         per = controls meaning of nve
c             = 0 for nve points per edge
c             = 1 for nve points per 2 pi on each edge,
c                 rounded up, so there is at least 1 point per edge;
c                 excess points are set to zero.
c         nve = desired number of points per edge
c             = 1 for vertices only
c             = 2 for vertices + midpoints of edges
c             >= 2 for vertices + (nve-1) non-vertex points per edge.
c         vcirc = 1 to return vertices and edge points also for bounding
c                 circles which do not intersect any other circle;
c               = 0 otherwise.
c Output: ve(3,ive,iv),ive=1,nve,iv=1,nv are points uniformly spaced
c                   along the edges of the polygon,
c                   starting with vertices in ve(3,1,iv).
c         angle(iv),i=1,nv are azimuthal lengths of edges of polygon.
c         ipv(iv),iv=1,nv is the circle number of the iv'th edge,
c                  that is,points ve(*,ive,iv) lie on circle number ipv.
c         gp(i),i=1,np = which group circle i belongs to:
c                circle i and circle j are friends if they intersect,
c                and friends of friends are friends.
c         ev(i),i=1,nev = end index of each connected sequence of edges;
c                 the number of non-zero elements of ev is the number
c                 nev of connected boundaries of the polygon,
c                 and the last nev0 of these are non-intersecting.
c         nv = number of points;
c              if this exceeds nvmax, then you should call gvert again
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
C     print *,'--------------------'
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
c        initialise iwk to inadmissible values
      do iv=1,nvmax
        iwk(iv,2)=-1
      enddo
      do iv=1,nvmax
        iwk(iv,3)=-1
      enddo
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
c........construct cartesian axes with z-axis along rp(i)
        call gaxisi(rp(1,i),xi,yi)
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
c        want vertices and midpoints for circles without intersections
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
              dph=php-phm
              nv=nv+1
              if (nv.le.nvmax) then
c        edge index of this vertex
                ipv(nv)=i
c        azimuthal length of edge
                angle(nv)=dph
c        edge points
                zv=1._10-cmi
                if (per.eq.0) then
                  mve=nve
                else
                  amve=nve*(dph/TWOPI)
                  mve=amve
                  if (dble(mve).lt.amve) mve=mve+1
                endif
                do ive=1,mve
                  if (cm(i).ge.0._10) then
                    ph=(phm*(mve-ive+1)+php*(ive-1))/dble(mve)
                  else
                    ph=(php*(mve-ive+1)+phm*(ive-1))/dble(mve)
                  endif
                  xv=si*cos(ph)
                  yv=si*sin(ph)
                  ve(1,ive,nv)=zv*rp(1,i)+xv*xi(1)+yv*yi(1)
                  ve(2,ive,nv)=zv*rp(2,i)+xv*xi(2)+yv*yi(2)
                  ve(3,ive,nv)=zv*rp(3,i)+xv*xi(3)+yv*yi(3)
                enddo
                do ive=mve+1,nve
                  ve(1,ive,nv)=0._10
                  ve(2,ive,nv)=0._10
                  ve(3,ive,nv)=0._10
                enddo
              endif
c        record endpoints of edge
              call gvtrail(scmi,np,ii,km,kp,iwk,iwk(1,2),nvmax,nv,
     *          ik,ikchk)
            enddo
            nev0=nev0+1
          endif
C         print *,'at',i,': full circle'
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
c             print *,'*** warning from gvert: near multiple intersectio
c    *n at',i,': edge',km,kp,' dph=',real(dph)
c             warn=.true.
c           endif
C           print *,'at',i,': edge',km,kp,
C    *        ' (',jm,' in',jml,jmu,',',jp,' in',jpl,jpu,' of',ni,')',
C    *        ' ph=',real(ph),' dph=',real(dph)
            nv=nv+1
            if (nv.le.nvmax) then
c        edge index of this vertex
              ipv(nv)=i
c        azimuthal length of edge
              angle(nv)=dph
c        edge points
              zv=1._10-cmi
              if (per.eq.0) then
                mve=nve
              else
                amve=nve*(dph/TWOPI)
                mve=amve
                if (dble(mve).lt.amve) mve=mve+1
              endif
              do ive=1,mve
                if (cm(i).ge.0._10) then
                  ph=(phm*(mve-ive+1)+php*(ive-1))/dble(mve)
                else
                  ph=(php*(mve-ive+1)+phm*(ive-1))/dble(mve)
                endif
                xv=si*cos(ph)
                yv=si*sin(ph)
                ve(1,ive,nv)=zv*rp(1,i)+xv*xi(1)+yv*yi(1)
                ve(2,ive,nv)=zv*rp(2,i)+xv*xi(2)+yv*yi(2)
                ve(3,ive,nv)=zv*rp(3,i)+xv*xi(3)+yv*yi(3)
              enddo
              do ive=mve+1,nve
                ve(1,ive,nv)=0._10
                ve(2,ive,nv)=0._10
                ve(3,ive,nv)=0._10
              enddo
            endif
c        record endpoints of edge
            call gvtrail(scmi,np,i,km,kp,iwk,iwk(1,2),nvmax,nv,ik,ikchk)
c        do another segment
          goto 220
        endif
  280 continue
c--------check on whether ik endpoints matched ki endpoints
      if (ikchk.ne.0._10) then
C       warn=.true.
C       print *,'*** from gvert: at tol =',tol,
C    *    ', ikchk=',ikchk,' should be 0'
c        retry with modified tolerance
        call gtol(tol,tolin)
        goto 100
      elseif (tol.gt.0._10) then
C       print *,'... from gvert: success at tol =',tol
      endif
c--------order vertices right-handedly about polygon
      if (nv.gt.0.and.nv.le.nvmax) then
c........which group of friends each circle belongs to
        call ggp(np,gp)
c........find order of vertices around polygon
        call gvord(np,nv,ipv,gp,ev,nev,iwk,iwk(1,2),iwk(1,3),iwk(1,4))
c........reorder vertices
        call vpermi(ipv,nv,iwk(1,2),iwk(1,3))
        call vpermdd(ve,3*nve,nv,iwk(1,2),wk)
        call vpermd(angle,nv,iwk(1,2),wk)
      endif
c--------done
  410 continue
C     if (warn) then
C       write (*,'(a3,a20,4a24)')
C    *    ' ','x','y','z','r','1-c'
C       write (*,'(i3,5g24.16)')
C    *    (j,(rp(i,j),i=1,3),sqrt(rp(1,j)**2+rp(2,j)**2+rp(3,j)**2),
C    *    cm(j),j=1,np)
C     endif
c     if (nv.gt.nvmax) then
c       print *,'*** from gvert: number of vertices =',nv,
c    *    ' exceeds maximum',nvmax
c     endif
      return
c
  420 print *,'*** from gvert: total failure at tol = ',tol
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
