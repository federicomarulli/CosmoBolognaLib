c-----------------------------------------------------------------------
c © A J S Hamilton 2001
c-----------------------------------------------------------------------
c *
c * Inline code to find the first nb ordered elements of an array.
c * The order is defined by the `function' order(i,j).
c * It's inlined to allow greater flexibility in defining the
c * order() function.  For example, the same code works for both
c * real and integer arrays.
c *
c * Returns index iord of these nb elements, ordered so that iord(1)
c * corresponds to the element having the most extreme value.
c * If nb .gt. na, then the last nb-na elements of iord are undefined.
c * Uses heapsort algorithm.
c *
      n=min(na,nb)
      if (n.le.0) return
c        heap first n elements, so smallest element is at top of heap
      do 210 ib=n/2+1,n
        iord(ib)=ib
  210 continue
      do 240 ia=n/2,1,-1
        i=ia
  230   ib=i*2
        if (ib.le.n) then
          if (ib.lt.n) then
            if (order(iord(ib),iord(ib+1))) ib=ib+1
          endif
          if (order(ia,iord(ib))) then
            iord(i)=iord(ib)
            i=ib
            goto 230
          endif
        endif
        iord(i)=ia
  240 continue
c        now compare rest of elements of array to heap
      do 340 ia=n+1,na
c        if new element is greater than smallest, sift it into heap
        i=1
        if (order(ia,iord(i))) then
  330     ib=i*2
          if (ib.le.n) then
            if (ib.lt.n) then
              if (order(iord(ib),iord(ib+1))) ib=ib+1
            endif
            if (order(ia,iord(ib))) then
              iord(i)=iord(ib)
              i=ib
              goto 330
            endif
          endif
          iord(i)=ia
        endif
  340 continue
c        unheap iord so largest element is at top
      do 390 ja=2,n
        ia=2+n-ja
        it=iord(ia)
        i=1
        iord(ia)=iord(i)
  380   ib=i*2
        if (ib.le.ia-1) then
          if (ib.lt.ia-1) then
            if (order(iord(ib),iord(ib+1))) ib=ib+1
          endif
          if (order(it,iord(ib))) then
            iord(i)=iord(ib)
            i=ib
            goto 380
          endif
        endif
        iord(i)=it
  390 continue
c
