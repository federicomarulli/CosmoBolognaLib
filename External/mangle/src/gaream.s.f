c-----------------------------------------------------------------------
c © A J S Hamilton 2001
c-----------------------------------------------------------------------
      subroutine gaream(area,areat,rp,cm,np,tol,verb,npg,npp,
     *  cmimin,cmimax,phi,iord,ldegen)
      integer np,verb,npg,npp,iord(np)
      logical ldegen
      real*10 area,areat,rp(3,np),cm(np),tol,cmimin,cmimax,phi(2,np)
c
c        parameters
      include 'pi.par'
      real*10 TWOPI
      parameter (TWOPI=2._10*PI)
c        intrinsics
      intrinsic abs
c *
c * Same as garea, but speed up matters by checking first whether
c * cm(npg) or cm(npp) lie outside |cmimin| and |cmimax|,
c * giving simpler result.
c * cmimin and cmimax are gotten from prior call to gcmlim.
c *
c * Note this subroutine assumes, without checking, that:
c * if np.eq.npg,
c *    then cm(npg).ge.0;
c * elseif np.eq.npp,
c *    then cm(npg).ge.0 & cm(npp).lt.0 & rp(.,npg)=rp(.,npp).
c *
c        sphere
      if (np.eq.npg) then
        if (cm(npg).le.abs(cmimin)) then
c        region excludes sphere
          if (cmimin.ge.0._10) then
            area=0._10
c        region encloses sphere
          elseif (cmimin.lt.0._10) then
            area=TWOPI*cm(npg)
          endif
        elseif (cm(npg).ge.abs(cmimax)) then
c        sphere encloses region
          if (cmimax.ge.0._10) then
            area=areat
c        sphere and region enclose each other
          elseif (cmimax.lt.0._10) then
            area=areat-TWOPI*(2._10-cm(npg))
          endif
c        sphere intersects boundary of region
        else
          call garea(area,rp,cm,np,tol,verb,phi,iord,ldegen)
        endif
c        annulus
      elseif (np.eq.npp) then
c        region is null
        if (cm(npg).le.-cm(npp)) then
          area=0._10
        elseif (cm(npg).le.abs(cmimin)) then
c        region excludes annulus
          if (cmimin.ge.0._10) then
            area=0._10
c        region encloses annulus
          elseif (cmimin.lt.0._10) then
            area=TWOPI*(cm(npg)+cm(npp))
          endif
        elseif (-cm(npp).ge.abs(cmimax)) then
c        annulus encloses region
          if (cmimax.ge.0._10) then
            area=0._10
c        annulus and region enclose each other
          elseif (cmimax.lt.0._10) then
            area=TWOPI*(cm(npg)+cm(npp))
          endif
        elseif (cm(npg).ge.abs(cmimax)
     *    .and.-cm(npp).le.abs(cmimin)) then
          if (cmimin.ge.0._10) then
c        annulus contains region
            if (cmimax.ge.0._10) then
              area=areat
c        outer ring of annulus and region enclose each other
            elseif (cmimax.lt.0._10) then
              area=areat-TWOPI*(2._10-cm(npg))
            endif
          elseif (cmimin.lt.0._10) then
c        inner ring of annulus and region enclose each other
            if (cmimax.ge.0._10) then
              area=areat-TWOPI*(2._10+cm(npp))
c        annulus and region enclose each other
            elseif (cmimax.lt.0._10) then
              area=areat-TWOPI*(2._10-cm(npg)-cm(npp))
            endif
          endif
c        annulus intersects boundary of region
        else
          call garea(area,rp,cm,np,tol,verb,phi,iord,ldegen)
        endif
c        np .ne. npg or npp
      else
        call garea(area,rp,cm,np,tol,verb,phi,iord,ldegen)
      endif
      return
      end
c
