c-----------------------------------------------------------------------
c © A J S Hamilton 2001
c-----------------------------------------------------------------------
      logical function gptin(rp,cm,np,rpi)
      integer np
      real*10 rp(3,np),cm(np),rpi(3)
c
c        intrinsics
      intrinsic abs
c        externals
      integer gzeroar
c        local (automatic) variables
      integer j
      real*10 cmij,cmj
c *
c * Determine whether unit direction rpi lies within region bounded by
c *    1 - r.rp(j) <= cm(j)  (if cm(j).ge.0)
c *    1 - r.rp(j) > -cm(j)  (if cm(j).lt.0)
c * for j=1,np where rp(j) are unit directions.
c *
c  Input: rp(3,j),j=1,np
c         cm(j),j=1,np
c         np
c         rpi(3)
c Output: gptin = .true. if point lies within region
c                 .false. if outside.
c
      gptin=.false.
c        check for point outside because one circle is null
      if (gzeroar(cm,np).eq.0) goto 410
c        check each boundary
      do 140 j=1,np
c        null boundary means no constraint
        if (cm(j).ge.2._10) goto 140
        cmj=abs(cm(j))
c        1-cos of angle between point and rp(j) direction
        cmij=((rpi(1)-rp(1,j))**2+(rpi(2)-rp(2,j))**2
     *    +(rpi(3)-rp(3,j))**2)/2._10
c        check if point is outside rp(j) boundary
        if (cm(j).ge.0._10) then
          if (cmij.gt.cmj) goto 410
        elseif (cm(j).lt.0._10) then
          if (cmij.le.cmj) goto 410
        endif
  140 continue
c        point survived all assails
      gptin=.true.
c        done
  410 continue
      return
      end
c
