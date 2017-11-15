c-----------------------------------------------------------------------
c © A J S Hamilton 2001
c-----------------------------------------------------------------------
      subroutine wlm(w,lmax1,im,nw,ri,phi,qphi,zi,ci,si,ph,dph,v)
      integer lmax1,im,nw,qphi
      real*10 w(im,nw),ri,phi,zi,ci,si,ph,dph,v(lmax1)
c
c        parameters
      real*10 HALF
      parameter (HALF=1._10/2._10)
      include 'pi.par'
      real*10 TWOPI
      parameter (TWOPI=2._10*PI)

c largest non-overflowing 2^(DBL_MAX_EXP-1) for real*10 (from values.h)
      integer DBL_MAX_EXP
      parameter (DBL_MAX_EXP=1024)

c        intrinsics
      intrinsic abs,max
c        local variables
      integer l,lm,l1,m,mmax1,mmin1,mn,mn1,mq,m1,n,nmax1,n1
      integer ed,edm,edn,ee,eem,een,ez,ezn,iez
      real*10 al,al1,am,an,cm,cmp,cmphi,cnph,cp,
     *  d,dim,dm,dme,dn,dpe,dre,d0,d1,d2,e,ed0,em,en,e1,e2,fn,q,qq,q1,
     *  smphi,snph,t,z,zm,zn,zp
      real*10 OVFLOW,UNFLOW
c *
c * Calculates contribution to spherical transform w(lm)
c * by boundary segment of window function w=sum w(lm)*Y(l,m).
c * See AJSH Notes V2 210.
c *
c * Definition: z(l,m)*sqrt(2*l+1)*exp(i*m*phi)=Y(l,m) (m.ge.0).
c * d(l,m,n) is matrix for rotation about cone y-axis which takes
c * the north pole (=z-axis) to cone north pole.
c *
c * Underflow:
c * Is a pain!  It only matters when lmax is large,
c * and in most cases the underflow doesn't matter anway,
c * and it requires a pile of ugly coding.
c * But there are cases where it matters.
c *
c * Here's a quick argument for the Legendre polynomials.
c * z(n,n) is [sin(th)]^n up to a factor of order 1,
c * and can underflow if sin(th) is small and n is large.
c * z(l,n) is obtained by recursion from z(n,n),
c * and is `not small' for |sin(th)| > n/l
c * (the classically allowed regime),
c * so underflow causes z(l,n) to be underestimated when
c * l > n/|sin(th)| = n/[epsilon^(1/n)],
c * where epsilon is a number small enough to underflow.
c * The smallest l occurs when n = - ln(epsilon), whereat sin(th) = 1/e
c * and l = e n, with e = 2.718... being the exponential.
c * If epsilon = 2^-1023 = 1e-308,
c * then the smallest l occurs when n = 709, whereat l = 1928.
c *
c * Numerical experiment indicates that underflow sets in earlier
c * for the rotation matrix -- about l = 200, for sin(th) = 1/e.
c *
c  Input: lmax1 = lmax+1 where lmax is maximum desired l of transform.
c         im = 1 means compute real part of harmonics only;
c              2 means compute both real and imaginary parts.
c              Note a region has pure real harmonics if it has mirror
c              symmetry through x-z plane.
c         nw = [(lmax+1)*(lmax+2)]/2
c         ri,phi,zi = axis of segment cone in cylindrical coordinates.
c         qphi: use phi+qphi*pi/2 in place of phi,
c               for better precision near phi = integer*pi/2 .
c         ci,si = cos and sin of opening angle of segment cone.
c         ph = azimuthal angle of centre point of segment,
c              zero being such that rotation from z-axis to cone axis
c              is about cone's y-axis.
c         dph = azimuthal angle subtended by segment.
c Input/Output: w(i,lm) = spherical transform
c            w(i,lm), i=1,im, lm=l*(l+1)/2+m+1, l=0,lmax, m=0,l;
c            w(1,lm) is real part, w(2,lm) is imaginary part;
c            w(i,lm) is ADDED to input w(i,lm).
c         Note w(l,-m)=(-)^m*[Complex conjugate of w(l,m)], just as
c              Y(l,-m)=(-)^m*[Complex conjugate of Y(l,m)].
c Work array: v should be dimensioned at least lmax1
c
      if (lmax1.le.0) goto 300
      if (dph.eq.0._10) goto 300
c        largest non-overflowing power of 2 for real*10
      OVFLOW=2._10**(DBL_MAX_EXP-1)
      UNFLOW=1._10/OVFLOW
      z=OVFLOW*0._10
      if (z.ne.0._10.or.UNFLOW.eq.0._10) then
        print *,'*** from wlm: DBL_MAX_EXP on your machine appears to be
     * <',DBL_MAX_EXP
        print *,' please modify DBL_MAX_EXP in wlm.s.f in the mangle dir
     *ectory, and recompile;'
        print *,' DBL_MAX_EXP may be defined in float.h somewhere on you
     *r machine'
        stop
      endif
c        cm=(zi-1)/2; cp=(zi+1)/2
c        done this way to ensure accuracy at small ri
      cm=-((1._10-zi)**2+ri**2)/4._10
      cp=((1._10+zi)**2+ri**2)/4._10
      nmax1=lmax1
      if (dph.eq.TWOPI.or.dph.eq.-TWOPI) nmax1=1
      if (si.eq.0._10) nmax1=1
      mmin1=1
      mmax1=lmax1
      do 280 n1=1,nmax1
        n=n1-1
        an=n
c--------calculate spherical harmonics v(l,n) in cone frame
        do l1=n1,lmax1
          al1=l1
          l=l1-1
          al=l
c........l = n = 0
          if (n.eq.0) then
            if (l.eq.0) then
c        zn=z(0,0)
              zn=1._10/sqrt(2._10*TWOPI)
              ezn=0
              ez=ezn
              fn=dph
c        z=z(0,0); zp=z(1,0)
              z=zn
              zp=ci*z
c        v(0,0)
              v(l1)=-fn*zp
c........l > 0, n = 0
            elseif (l.gt.0) then
c        zm=z(l-1,0); z=z(l,0); zp=z(l+1,0)
              qq=al1+al
              zm=z
              z=zp
              zp=(qq*ci*z-al*zm)/al1
c        v(l,0)
              v(l1)=fn/sqrt(qq)*(zm-zp)
            endif
c........l = n > 0
          elseif (l.eq.n) then
c        zn=z(l,l)
            t=-sqrt((an-HALF)/an)*si
            z=t*zn
c        zn may underflow
            if (abs(z).le.UNFLOW) then
              zn=zn*OVFLOW
              z=t*zn
              ezn=ezn+1
            endif
            ez=ezn
            zn=z
            fn=2._10*sin(an*dph/2._10)/an
c        z=z(l,l); zp=z(l+1,l)
            q=sqrt(al1+al)
            zp=q*ci*z
c        v=v(l,l)
            v(l1)=-fn/al1*zp
c........l > n > 0
          elseif (l.gt.n) then
c        zm=z(l-1,n); z=z(l,n); zp=z(l+1,n)
            qq=al1+al
            q1=q
            q=sqrt((al1+an)*(al1-an))
            zm=z
            z=zp
            zp=(qq*ci*z-q1*zm)/q
            if (ez.gt.0) then
c        recover from underflow
              if (abs(zp).ge.1._10
     *          .and.(abs(z).ge.1._10.or.z.eq.0._10)
     *          .and.(abs(zm).ge.1._10.or.zm.eq.0._10)) then
                ez=ez-1
                zm=zm*UNFLOW
                z=z*UNFLOW
                zp=zp*UNFLOW
              endif
            endif
c        v=v(l,n)
            v(l1)=fn/sqrt(qq)*(q1/al*zm-q/al1*zp)
          endif
c........restore correct scaling
          if (ez.gt.0) then
            do iez=1,ez
              v(l1)=v(l1)*UNFLOW
            enddo
          endif
        enddo
c........skip rotation if harmonics in cone frame are all zero
        do l1=lmax1,n1,-1
          if (v(l1).ne.0._10) goto 220
        enddo
        goto 280
  220   continue
c--------matrix d(l,m,n) to rotate v(l,n) about cone y-axis
        cnph=cos(an*ph)
        snph=sin(an*ph)
c        cone axis is parallel to desired axis
        if (ri.eq.0._10) then
          mmin1=n1
          mmax1=n1
        endif
        do m1=mmin1,mmax1
          m=m1-1
          am=m
          cmphi=cos(am*phi)
          smphi=sin(am*phi)
          if (qphi.ne.0) then
            mq=m*qphi
            if (mod(mq/2,2).ne.0) then
              cmphi=-cmphi
              smphi=-smphi
            endif
            if (mod(mq,2).eq.1) then
              cmp=cmphi
              cmphi=-smphi
              smphi=cmp
            elseif (mod(mq,2).eq.-1) then
              cmp=cmphi
              cmphi=smphi
              smphi=-cmp
            endif
          endif
c........m = n = 0
          if (m.eq.0) then
            if (n.eq.0) then
c        d0=d(0,0,0)
              d0=1._10
              ed0=0
c        initialize dn, en
              dn=d0
              en=d0
              edn=ed0
              een=ed0
c........m = 0, n > 0
            elseif (n.gt.0) then
c        d0=d(n,0,n)=(-)^n*d(n,0,-n)
              t=sqrt((an-HALF)/an)*ri
              d=d0*t
c        d0 may underflow
              if (abs(d).le.UNFLOW) then
                d0=d0*OVFLOW
                d=d0*t
                ed0=ed0+1
              endif
              d0=d
            endif
c        initialize dm, em
            dm=d0
            em=d0
            edm=ed0
            eem=ed0
c........0 < m < n
          elseif (m.lt.n) then
c        dm=d(n,m,n); em=(-)^n*d(n,m,-n)
            t=sqrt((an-am+1._10)/(an+am))*2._10/ri
c        dm
            d=dm*cp*t
c        dm may underflow ...
            if (abs(d).le.UNFLOW.and.dm.ne.0._10) then
              dm=dm*OVFLOW
              d=dm*cp*t
              edm=edm+1
c        ... or may recover from underflow
            elseif (edm.gt.0) then
              if (abs(d).ge.1._10) then
                d=d*UNFLOW
                edm=edm-1
              endif
            endif
            dm=d
c        em
            e=em*cm*t
c        em may underflow ...
            if (abs(e).le.UNFLOW.and.em.ne.0._10) then
              em=em*OVFLOW
              e=em*cm*t
              eem=eem+1
c        ... or may recover from underflow
            elseif (eem.gt.0) then
              if (abs(e).ge.1._10) then
                e=e*UNFLOW
                eem=eem-1
              endif
            endif
            em=e
c........0 < m = n
          elseif (m.eq.n) then
c        dn=d(n,n,n); en=(-)^n*d(n,n,-n)
c        dn
            d=cp*dn
c        dn may underflow
            if (abs(d).le.UNFLOW.and.dn.ne.0._10) then
              dn=dn*OVFLOW
              d=cp*dn
              edn=edn+1
            endif
            dn=d
c        en
            e=cm*en
c        en may underflow
            if (abs(e).le.UNFLOW.and.en.ne.0._10) then
              en=en*OVFLOW
              e=cm*en
              een=een+1
            endif
            en=e
c        initialize dm, em
            dm=dn
            em=en
            edm=edn
            eem=een
c........m > n
          elseif (m.gt.n) then
c        dm=d(m,m,n); em=(-)^n*d(m,m,-n)
            t=-sqrt(am*(am-HALF)/((am+an)*(am-an)))*ri
c        dm
            d=t*dm
c        dm may underflow ...
            if (abs(d).le.UNFLOW.and.dm.ne.0._10) then
              dm=dm*OVFLOW
              d=t*dm
              edm=edm+1
c        ... or may recover from underflow
            elseif (edm.gt.0) then
              if (abs(d).ge.1._10) then
                d=d*UNFLOW
                edm=edm-1
              endif
            endif
            dm=d
c        em
            if (n.gt.0) then
              e=t*em
c        em may underflow ...
              if (abs(e).le.UNFLOW.and.em.ne.0._10) then
                em=em*OVFLOW
                e=t*em
                eem=eem+1
c        ... or may recover from underflow
              elseif (eem.gt.0) then
                if (abs(e).ge.1._10) then
                  e=e*UNFLOW
                  eem=eem-1
                endif
              endif
              em=e
            endif
          endif
c........l >= m, n
          mn=max(m,n)
          mn1=mn+1
          do 240 l1=mn1,lmax1
            l=l1-1
            al=l
c        d=d(l,m,n); e=(-)^n*d(l,m,-n)
            if (l.eq.mn) then
              lm=(l*l1)/2+m1
              q=0._10
              d1=0._10
              d=dm
              ed=edm
              if (n.eq.0) then
                e=0._10
                ee=0
              elseif (n.gt.0) then
                e1=0._10
                e=em
                ee=eem
              endif
            elseif (l.gt.mn) then
              lm=lm+l
              qq=2._10*al-1._10
              q1=q
              q=sqrt((al+am)*(al-am)*(al+an)*(al-an))/al
              t=0._10
              if (m*n.gt.0) t=am*an/((al-1._10)*al)
              d2=d1
              d1=d
              d=(qq*(zi-t)*d1-q1*d2)/q
c        recover d from underflow
              if (ed.gt.0) then
                if (abs(d).ge.1._10
     *            .and.(abs(d1).ge.1._10.or.e1_10.eq.0._10)
     *            .and.(abs(d2).ge.1._10.or.e2_10.eq.0._10)) then
                  d2=d2*UNFLOW
                  d1=d1*UNFLOW
                  d=d*UNFLOW
                  ed=ed-1
                endif
              endif
              if (n.gt.0) then
                e2=e1
                e1=e
                e=(qq*(zi+t)*e1-q1*e2)/q
c        recover e from underflow
                if (ee.gt.0) then
                  if (abs(e).ge.1._10
     *              .and.(abs(e1).ge.1._10.or.e1.eq.0._10)
     *              .and.(abs(e2).ge.1._10.or.e2.eq.0._10)) then
                    e2=e2*UNFLOW
                    e1=e1*UNFLOW
                    e=e*UNFLOW
                    ee=ee-1
                  endif
                endif
              endif
            endif
c--------rotate: w(l,m) = d(l,m,0)*v(l,0) +
c            sum (n=1,l) [d(l,m,n)*v(l,n) + d(l,m,-n)*v(l,-n)]
c           if ((l1.eq.lmax1.or.l1.eq.n1)
c    *        .and.(m1.eq.n1.or.m1.eq.1)) then
c             write (*,'("l m n = ",3i4," d(l,m,n) = ",
c    *          2g13.5,i3,"   (-)^n d(l,m,-n) = ",2g13.5,i3)')
c    *          l1-1,m1-1,n1-1,d,d*UNFLOW**ed,ed,e,e*UNFLOW**ee,ee
c           endif
            if (ed.eq.0) then
              if (ee.eq.0) then
                dpe=d+e
                dme=d-e
              else
                dpe=d
                dme=d
              endif
            else
              if (ee.eq.0) then
                dpe=e
                dme=-e
              else
                goto 240
              endif
            endif
            dre=cnph*dpe*v(l1)
            dim=snph*dme*v(l1)
            w(1,lm)=w(1,lm)+cmphi*dre-smphi*dim
            if (im.eq.2) w(2,lm)=w(2,lm)-cmphi*dim-smphi*dre
  240     continue
        enddo
  280 continue
c--------done
  300 continue
      return
      end
c
