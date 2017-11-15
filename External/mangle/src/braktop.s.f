c-----------------------------------------------------------------------
c © A J S Hamilton 2001
c-----------------------------------------------------------------------
      subroutine braktop(aa,ia,a,n,l)
      integer ia,n,l
      real*10 aa,a(n)
      integer idir,istep
c *
c * Bracket aa in table a ordered in decreasing order,
c * i.e. find ia such that
c *    a(ia) >= aa > a(ia+1)   (l.eq.0)
c * or a(ia) > aa >= a(ia+1)   (l.eq.1).
c * Gives ia=0 if
c *    aa > a(1)   (l.eq.0)
c * or aa >= a(1)  (l.eq.1).
c * Gives ia=n if
c *    aa <= a(1)  (l.eq.0)
c * or aa < a(1)   (l.eq.1).
c *
c  Input: aa = value to bracket
c         a = ordered array of length n
c         n = length of a array
c         l = 0 or 1 to choose equality on lower or upper limit.
c Input/Output: ia on input = element of a to start search from;
c                  on output = as above.
c
      if (ia.lt.1) ia=1
      if (ia.gt.n) ia=n
c        idir is which way to search, up or down
      if (aa.lt.a(ia)) then
        idir=1
      elseif (aa.gt.a(ia)) then
        idir=-1
      elseif (aa.eq.a(ia)) then
        ia=ia-l
        goto 200
      endif
c        istep is how far to leap from present position
      istep=1
c        leap
  120 ia=ia+idir*istep
c        keep doubling leap till you've straddled desired place
      if (ia.lt.1.or.ia.gt.n) then
        continue
      elseif (idir*aa.lt.idir*a(ia)) then
        istep=istep*2
        goto 120
      elseif (aa.eq.a(ia)) then
        ia=ia-l
        goto 200
      endif
      idir=-idir
c        binary chop homes in on desired place in table
  140 if (istep.gt.1) then
        istep=istep/2
        ia=ia+idir*istep
        if (ia.lt.1) then
          idir=1
        elseif (ia.gt.n) then
          idir=-1
        elseif (aa.lt.a(ia)) then
          idir=1
        elseif (aa.gt.a(ia)) then
          idir=-1
        elseif (aa.eq.a(ia)) then
          ia=ia-l
          goto 200
        endif
        goto 140
      endif
      if (idir.eq.-1) ia=ia-1
  200 continue
      return
      end
c
c-----------------------------------------------------------------------
      subroutine brakbot(aa,ia,a,n,l)
      integer ia,n,l
      real*10 aa,a(n)
      integer idir,istep
c *
c * Bracket aa in table a ordered in increasing order,
c * i.e. find ia such that
c *    a(ia) <= aa < a(ia+1)   (l.eq.0)
c * or a(ia) < aa <= a(ia+1)   (l.eq.1).
c * Gives ia=0 if
c *    aa < a(1)   (l.eq.0)
c * or aa <= a(1)  (l.eq.1).
c * Gives ia=n if
c *    aa <= a(1)  (l.eq.0)
c * or aa < a(1)   (l.eq.1).
c *
c  Input: aa = value to bracket
c         a = ordered array of length n
c         n = length of a array
c         l = 0 or 1 to choose equality on lower or upper limit.
c Input/Output: ia on input = element of a to start search from;
c                  on output = as above.
c
      if (ia.lt.1) ia=1
      if (ia.gt.n) ia=n
c        idir is which way to search, up or down
      if (aa.gt.a(ia)) then
        idir=1
      elseif (aa.lt.a(ia)) then
        idir=-1
      elseif (aa.eq.a(ia)) then
        ia=ia-l
        goto 200
      endif
c        istep is how far to leap from present position
      istep=1
c        leap
  120 ia=ia+idir*istep
c        keep doubling leap till you've straddled desired place
      if (ia.lt.1.or.ia.gt.n) then
        continue
      elseif (idir*aa.gt.idir*a(ia)) then
        istep=istep*2
        goto 120
      elseif (aa.eq.a(ia)) then
        ia=ia-l
        goto 200
      endif
      idir=-idir
c        binary chop homes in on desired place in table
  140 if (istep.gt.1) then
        istep=istep/2
        ia=ia+idir*istep
        if (ia.lt.1) then
          idir=1
        elseif (ia.gt.n) then
          idir=-1
        elseif (aa.gt.a(ia)) then
          idir=1
        elseif (aa.lt.a(ia)) then
          idir=-1
        elseif (aa.eq.a(ia)) then
          ia=ia-l
          goto 200
        endif
        goto 140
      endif
      if (idir.eq.-1) ia=ia-1
  200 continue
      return
      end
c
c-----------------------------------------------------------------------
      subroutine braktpa(aa,ia,a,n,l)
      integer ia,n,l
      real*10 aa,a(n)
c        intrinsics
      intrinsic abs
c        local (automatic) variables
      integer idir,istep
c *
c * Bracket aa in table a ordered in decreasing order of absolute value,
c * i.e. find ia such that
c *    |a(ia)| >= aa > |a(ia+1)|   (l.eq.0)
c * or |a(ia)| > aa >= |a(ia+1)|   (l.eq.1).
c * Gives ia=0 if
c *    aa > |a(1)|   (l.eq.0)
c * or aa >= |a(1)|  (l.eq.1).
c * Gives ia=n if
c *    aa <= |a(1)|  (l.eq.0)
c * or aa < |a(1)|   (l.eq.1).
c *
c  Input: aa = value to bracket
c         a = ordered array of length n
c         n = length of a array
c         l = 0 or 1 to choose equality on lower or upper limit.
c Input/Output: ia on input = element of a to start search from;
c                  on output = as above.
c
      if (ia.lt.1) ia=1
      if (ia.gt.n) ia=n
c        idir is which way to search, up or down
      if (aa.lt.abs(a(ia))) then
        idir=1
      elseif (aa.gt.abs(a(ia))) then
        idir=-1
      elseif (aa.eq.abs(a(ia))) then
        ia=ia-l
        goto 200
      endif
c        istep is how far to leap from present position
      istep=1
c        leap
  120 ia=ia+idir*istep
c        keep doubling leap till you've straddled desired place
      if (ia.lt.1.or.ia.gt.n) then
        continue
      elseif (idir*aa.lt.idir*abs(a(ia))) then
        istep=istep*2
        goto 120
      elseif (aa.eq.abs(a(ia))) then
        ia=ia-l
        goto 200
      endif
      idir=-idir
c        binary chop homes in on desired place in table
  140 if (istep.gt.1) then
        istep=istep/2
        ia=ia+idir*istep
        if (ia.lt.1) then
          idir=1
        elseif (ia.gt.n) then
          idir=-1
        elseif (aa.lt.abs(a(ia))) then
          idir=1
        elseif (aa.gt.abs(a(ia))) then
          idir=-1
        elseif (aa.eq.abs(a(ia))) then
          ia=ia-l
          goto 200
        endif
        goto 140
      endif
      if (idir.eq.-1) ia=ia-1
  200 continue
      return
      end
c
c-----------------------------------------------------------------------
      subroutine brakbta(aa,ia,a,n,l)
      integer ia,n,l
      real*10 aa,a(n)
c        intrinsics
      intrinsic abs
c        local (automatic) variables
      integer idir,istep
c *
c * Bracket aa in table a ordered in increasing order of absolute value,
c * i.e. find ia such that
c *    |a(ia)| <= aa < |a(ia+1)|   (l.eq.0)
c * or |a(ia)| < aa <= |a(ia+1)|   (l.eq.1).
c * Gives ia=0 if
c *    aa < |a(1)|   (l.eq.0)
c * or aa <= |a(1)|  (l.eq.1).
c * Gives ia=n if
c *    aa <= |a(1)|  (l.eq.0)
c * or aa < |a(1)|   (l.eq.1).
c *
c  Input: aa = value to bracket
c         a = ordered array of length n
c         n = length of a array
c         l = 0 or 1 to choose equality on lower or upper limit.
c Input/Output: ia on input = element of a to start search from;
c                  on output = as above.
c
      if (ia.lt.1) ia=1
      if (ia.gt.n) ia=n
c        idir is which way to search, up or down
      if (aa.gt.abs(a(ia))) then
        idir=1
      elseif (aa.lt.abs(a(ia))) then
        idir=-1
      elseif (aa.eq.abs(a(ia))) then
        ia=ia-l
        goto 200
      endif
c        istep is how far to leap from present position
      istep=1
c        leap
  120 ia=ia+idir*istep
c        keep doubling leap till you've straddled desired place
      if (ia.lt.1.or.ia.gt.n) then
        continue
      elseif (idir*aa.gt.idir*abs(a(ia))) then
        istep=istep*2
        goto 120
      elseif (aa.eq.abs(a(ia))) then
        ia=ia-l
        goto 200
      endif
      idir=-idir
c        binary chop homes in on desired place in table
  140 if (istep.gt.1) then
        istep=istep/2
        ia=ia+idir*istep
        if (ia.lt.1) then
          idir=1
        elseif (ia.gt.n) then
          idir=-1
        elseif (aa.gt.abs(a(ia))) then
          idir=1
        elseif (aa.lt.abs(a(ia))) then
          idir=-1
        elseif (aa.eq.abs(a(ia))) then
          ia=ia-l
          goto 200
        endif
        goto 140
      endif
      if (idir.eq.-1) ia=ia-1
  200 continue
      return
      end
c
