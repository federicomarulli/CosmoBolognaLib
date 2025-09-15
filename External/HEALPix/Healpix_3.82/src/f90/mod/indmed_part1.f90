! to be included in indmed.f90 (part 1/2)
      NMED = (NDAT+1) / 2
      IWRKT(1:ndat) = IDATT(:)
!
!  If the number of values is small, then use insertion sort
!
     If (NDAT < 35) Then
!
!  Bring minimum to first location to save test in decreasing loop
!
        IDCR = NDAT
        If (XDATT (IWRKT (1)) < XDATT (IWRKT (IDCR))) Then
           IWRK = IWRKT (1)
        Else
           IWRK = IWRKT (IDCR)
           IWRKT (IDCR) = IWRKT (1)
        Endif
        XWRK = XDATT (IWRK)
        Do ITMP = 1, NDAT - 2
           IDCR = IDCR - 1
           IWRK1 = IWRKT (IDCR)
           XWRK1 = XDATT (IWRK1)
           If (XWRK1 < XWRK) Then
              IWRKT (IDCR) = IWRK
              XWRK = XWRK1
              IWRK = IWRK1
           Endif
        End Do
        IWRKT (1) = IWRK
!
! Sort the first half, until we have NMED sorted values
!
        Do ICRS = 3, NMED
           XWRK = XDATT (IWRKT (ICRS))
           IWRK = IWRKT (ICRS)
           IDCR = ICRS - 1
           Do
              If (XWRK >= XDATT (IWRKT(IDCR))) Exit
              IWRKT (IDCR+1) = IWRKT (IDCR)
              IDCR = IDCR - 1
           End Do
           IWRKT (IDCR+1) = IWRK
        End Do
!
!  Insert any value less than the current median in the first half
!
        XWRK1 = XDATT (IWRKT (NMED))
        Do ICRS = NMED+1, NDAT
           XWRK = XDATT (IWRKT (ICRS))
           IWRK = IWRKT (ICRS)
           If (XWRK < XWRK1) Then
              IDCR = NMED - 1
              Do
                 If (XWRK >= XDATT (IWRKT(IDCR))) Exit
                 IWRKT (IDCR+1) = IWRKT (IDCR)
                 IDCR = IDCR - 1
              End Do
              IWRKT (IDCR+1) = IWRK
              XWRK1 = XDATT (IWRKT (NMED))
           End If
        End Do
        ires_med = IWRKT (NMED)
        Return
     End If
!
!  Make sorted subsets of 7 elements
!  This is done by a variant of insertion sort where a first
!  pass is used to bring the smallest element to the first position
!  decreasing disorder at the same time, so that we may remove
!  remove the loop test in the insertion loop.
!
     IMAX = 1
     IMIN = 1
     XMAX = XDATT (IWRKT(IMAX))
     XMIN = XDATT (IWRKT(IMIN))
     DO IDEB = 1, NDAT-6, 7
        IDCR = IDEB + 6
        If (XDATT (IWRKT(IDEB)) < XDATT (IWRKT(IDCR))) Then
           IWRK = IWRKT(IDEB)
        Else
           IWRK = IWRKT (IDCR)
           IWRKT (IDCR) = IWRKT(IDEB)
        Endif
        XWRK = XDATT (IWRK)
        Do ITMP = 1, 5
           IDCR = IDCR - 1
           IWRK1 = IWRKT (IDCR)
           XWRK1 = XDATT (IWRK1)
           If (XWRK1 < XWRK) Then
              IWRKT (IDCR) = IWRK
              IWRK = IWRK1
              XWRK = XWRK1
           Endif
        End Do
        IWRKT (IDEB) = IWRK
        If (XWRK < XMIN) Then
           IMIN = IWRK
           XMIN = XWRK
        End If
        Do ICRS = IDEB+1, IDEB+5
           IWRK = IWRKT (ICRS+1)
           XWRK = XDATT (IWRK)
           IDON = IWRKT(ICRS)
           If (XWRK < XDATT(IDON)) Then
              IWRKT (ICRS+1) = IDON
              IDCR = ICRS
              IWRK1 = IWRKT (IDCR-1)
              XWRK1 = XDATT (IWRK1)
              Do
                 If (XWRK >= XWRK1) Exit
                 IWRKT (IDCR) = IWRK1
                 IDCR = IDCR - 1
                 IWRK1 = IWRKT (IDCR-1)
                 XWRK1 = XDATT (IWRK1)
              End Do
              IWRKT (IDCR) = IWRK
           EndIf
        End Do
        If (XWRK > XMAX) Then
           IMAX = IWRK
           XMAX = XWRK
        End If
     End Do
!
!  Add-up alternatively MAX and MIN values to make the number of data
!  an exact multiple of 7.
!
     IDEB = 7 * (NDAT/7)
     NTRI = NDAT
     If (IDEB < NDAT) Then
!
        Do ICRS = IDEB+1, NDAT
           XWRK1 = XDATT (IWRKT (ICRS))
           IF (XWRK1 > XMAX) Then
              IMAX = IWRKT (ICRS)
              XMAX = XWRK1
           End If
           IF (XWRK1 < XMIN) Then
              IMIN = IWRKT (ICRS)
              XMIN = XWRK1
           End If
        End Do
        IWRK1 = IMAX
        Do ICRS = NDAT+1, IDEB+7
           IWRKT (ICRS) = IWRK1
           If (IWRK1 == IMAX) Then
              IWRK1 = IMIN
           Else
              NMED = NMED + 1
              IWRK1 = IMAX
           End If
        End Do
!
        Do ICRS = IDEB+2, IDEB+7
           IWRK = IWRKT (ICRS)
           XWRK = XDATT (IWRK)
           Do IDCR = ICRS - 1, IDEB+1, - 1
              If (XWRK >= XDATT (IWRKT(IDCR))) Exit
              IWRKT (IDCR+1) = IWRKT (IDCR)
           End Do
           IWRKT (IDCR+1) = IWRK
        End Do
!
        NTRI = IDEB+7
     End If
!
!  Make the set of the indices of median values of each sorted subset
!
     IDON1 = 0
     Do IDON = 1, NTRI, 7
        IDON1 = IDON1 + 1
        IMEDT (IDON1) = IWRKT (IDON + 3)
     End Do
!
!  Find XMED7, the median of the medians
!
