c-----------------------------------------------------------------------
c © A J S Hamilton 2001
c-----------------------------------------------------------------------
      subroutine findtop(a,na,iord,nb)
      integer na,nb,iord(nb)
      real*10 a(na)
c
c        local (automatic) variables
      integer i,ia,ib,it,n,ja
      logical order
c *
c * Find nb elements of real*10 a having the largest value.
c * Returns index iord of these elements, ordered so iord(1) corresponds
c * to element a(iord(1)) having largest value.
c * If nb .gt. na, last nb-na elements of iord are undefined.
c * Elements of a that are equal are left in their original order.
c *
      order(ia,ja)=a(ia).gt.a(ja).or.(a(ia).eq.a(ja).and.ia.lt.ja)
      include 'heapsort.inc'
      return
      end
c
c-----------------------------------------------------------------------
      subroutine findbot(a,na,iord,nb)
      integer na,nb,iord(nb)
      real*10 a(na)
c
c        local (automatic) variables
      integer i,ia,ib,it,n,ja
      logical order
c *
c * Find nb elements of real*10 a having the smallest value.
c * Returns index iord of these elements, ordered so iord(1) corresponds
c * to element a(iord(1)) having smallest value.
c * If nb .gt. na, last nb-na elements of iord are undefined.
c * Elements of a that are equal are left in their original order.
c *
      order(ia,ja)=a(ia).lt.a(ja).or.(a(ia).eq.a(ja).and.ia.lt.ja)
      include 'heapsort.inc'
      return
      end
c
c-----------------------------------------------------------------------
      subroutine findtpa(a,na,iord,nb)
      integer na,nb,iord(nb)
      real*10 a(na)
c
c        intrinsics
      intrinsic abs
c        local (automatic) variables
      integer i,ia,ib,it,n,ja
      logical order
c *
c * Find nb elements of real*10 a having the largest absolute value.
c * Returns index iord of these elements, ordered so iord(1) corresponds
c * to element a(iord(1)) having largest absolute value.
c * If nb .gt. na, last nb-na elements of iord are undefined.
c * Elements of a equal in abs value are left in their original order.
c *
      order(ia,ja)=abs(a(ia)).gt.abs(a(ja))
     *  .or.(abs(a(ia)).eq.abs(a(ja)).and.ia.lt.ja)
      include 'heapsort.inc'
      return
      end
c
c-----------------------------------------------------------------------
      subroutine findbta(a,na,iord,nb)
      integer na,nb,iord(nb)
      real*10 a(na)
c
c        intrinsics
      intrinsic abs
c        local (automatic) variables
      integer i,ia,ib,it,n,ja
      logical order
c *
c * Find nb elements of real*10 a having the smallest absolute value.
c * Returns index iord of these elements, ordered so iord(1) corresponds
c * to element a(iord(1)) having smallest absolute value.
c * If nb .gt. na, last nb-na elements of iord are undefined.
c * Elements of a equal in abs value are left in their original order.
c *
      order(ia,ja)=abs(a(ia)).lt.abs(a(ja))
     *  .or.(abs(a(ia)).eq.abs(a(ja)).and.ia.lt.ja)
      include 'heapsort.inc'
      return
      end
c
c-----------------------------------------------------------------------
      subroutine finitop(a,na,iord,nb)
      integer na,nb,iord(nb)
      integer a(na)
c
c        local (automatic) variables
      integer i,ia,ib,it,n,ja
      logical order
c *
c * Find nb elements of integer a having the largest value.
c * Returns index iord of these elements, ordered so iord(1) corresponds
c * to element a(iord(1)) having largest value.
c * If nb .gt. na, last nb-na elements of iord are undefined.
c * Elements of a that are equal are left in their original order.
c *
      order(ia,ja)=a(ia).gt.a(ja).or.(a(ia).eq.a(ja).and.ia.lt.ja)
      include 'heapsort.inc'
      return
      end
c
c-----------------------------------------------------------------------
      subroutine finibot(a,na,iord,nb)
      integer na,nb,iord(nb)
      integer a(na)
c
c        local (automatic) variables
      integer i,ia,ib,it,n,ja
      logical order
c *
c * Find nb elements of integer a having the smallest value.
c * Returns index iord of these elements, ordered so iord(1) corresponds
c * to element a(iord(1)) having smallest value.
c * If nb .gt. na, last nb-na elements of iord are undefined.
c * Elements of a that are equal are left in their original order.
c *
      order(ia,ja)=a(ia).lt.a(ja).or.(a(ia).eq.a(ja).and.ia.lt.ja)
      include 'heapsort.inc'
      return
      end
c
c-----------------------------------------------------------------------
      subroutine finitpa(a,na,iord,nb)
      integer na,nb,iord(nb)
      integer a(na)
c
c        intrinsics
      intrinsic abs
c        local (automatic) variables
      integer i,ia,ib,it,n,ja
      logical order
c *
c * Find nb elements of integer a having the largest absolute value.
c * Returns index iord of these elements, ordered so iord(1) corresponds
c * to element a(iord(1)) having largest absolute value.
c * If nb .gt. na, last nb-na elements of iord are undefined.
c * Elements of a equal in abs value are left in their original order.
c *
      order(ia,ja)=abs(a(ia)).gt.abs(a(ja))
     *  .or.(abs(a(ia)).eq.abs(a(ja)).and.ia.lt.ja)
      include 'heapsort.inc'
      return
      end
c
c-----------------------------------------------------------------------
      subroutine finibta(a,na,iord,nb)
      integer na,nb,iord(nb)
      integer a(na)
c
c        intrinsics
      intrinsic abs
c        local (automatic) variables
      integer i,ia,ib,it,n,ja
      logical order
c *
c * Find nb elements of integer a having the smallest absolute value.
c * Returns index iord of these elements, ordered so iord(1) corresponds
c * to element a(iord(1)) having smallest absolute value.
c * If nb .gt. na, last nb-na elements of iord are undefined.
c * Elements of a equal in abs value are left in their original order.
c *
      order(ia,ja)=abs(a(ia)).lt.abs(a(ja))
     *  .or.(abs(a(ia)).eq.abs(a(ja)).and.ia.lt.ja)
      include 'heapsort.inc'
      return
      end
c
