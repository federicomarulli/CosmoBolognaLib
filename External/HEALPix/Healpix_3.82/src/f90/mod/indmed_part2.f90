! to be included in indmed.f90 (part 2/2)
     XMED7 = XDATT (IMED7)
!
!  Count how many values are not higher than (and how many equal to) XMED7
!  This number is at least 4 * 1/2 * (N/7) : 4 values in each of the
!  subsets where the median is lower than the median of medians. For similar
!  reasons, we also have at least 2N/7 values not lower than XMED7. At the
!  same time, we find in each subset the index of the last value < XMED7,
!  and that of the first > XMED7. These indices will be used to restrict the
!  search for the median as the Kth element in the subset (> or <) where
!  we know it to be.
!
     IDON1 = 1
     NLEQ = 0
     NEQU = 0
     Do IDON = 1, NTRI, 7
        IMED = IDON+3
        If (XDATT (IWRKT (IMED)) > XMED7) Then
           IMED = IMED - 2
           If (XDATT (IWRKT (IMED)) > XMED7) Then
              IMED = IMED - 1
           Else If (XDATT (IWRKT (IMED)) < XMED7) Then
              IMED = IMED + 1
           Endif
        Else If (XDATT (IWRKT (IMED)) < XMED7) Then
           IMED = IMED + 2
           If (XDATT (IWRKT (IMED)) > XMED7) Then
              IMED = IMED - 1
           Else If (XDATT (IWRKT (IMED)) < XMED7) Then
              IMED = IMED + 1
           Endif
        Endif
        If (XDATT (IWRKT (IMED)) > XMED7) Then
           NLEQ = NLEQ + IMED - IDON
           IENDT (IDON1) = IMED - 1
           ISTRT (IDON1) = IMED
        Else If (XDATT (IWRKT (IMED)) < XMED7) Then
           NLEQ = NLEQ + IMED - IDON + 1
           IENDT (IDON1) = IMED
           ISTRT (IDON1) = IMED + 1
        Else                    !       If (XDATT (IWRKT (IMED)) == XMED7)
           NLEQ = NLEQ + IMED - IDON + 1
           NEQU = NEQU + 1
           IENDT (IDON1) = IMED - 1
           Do IMED1 = IMED - 1, IDON, -1
              If (XDATT (IWRKT (IMED1)) == XMED7) Then
                 NEQU = NEQU + 1
                 IENDT (IDON1) = IMED1 - 1
              Else
                 Exit
              End If
           End Do
           ISTRT (IDON1) = IMED + 1
           Do IMED1 = IMED + 1, IDON + 6
              If (XDATT (IWRKT (IMED1)) == XMED7) Then
                 NEQU = NEQU + 1
                 NLEQ = NLEQ + 1
                 ISTRT (IDON1) = IMED1 + 1
              Else
                 Exit
              End If
           End Do
        Endif
        IDON1 = IDON1 + 1
     End Do
!
!  Carry out a partial insertion sort to find the Kth smallest of the
!  large values, or the Kth largest of the small values, according to
!  what is needed.
!
!
     If (NLEQ - NEQU + 1 <= NMED) Then
        If (NLEQ < NMED) Then   !      Not enough low values
           IWRK1 = IMAX
           XWRK1 = XDATT (IWRK1)
           NORD = NMED - NLEQ
           IDON1 = 0
           ICRS1 = 1
           ICRS2 = 0
           IDCR = 0
           Do IDON = 1, NTRI, 7
              IDON1 = IDON1 + 1
              If (ICRS2 < NORD) Then
                 Do ICRS = ISTRT (IDON1), IDON + 6
                    If (XDATT (IWRKT (ICRS)) < XWRK1) Then
                       IWRK = IWRKT (ICRS)
                       XWRK = XDATT (IWRK)
                       Do IDCR = ICRS1 - 1, 1, - 1
                          If (XWRK >= XDATT (IWRKT (IDCR))) Exit
                          IWRKT  (IDCR+1) = IWRKT (IDCR)
                       End Do
                       IWRKT (IDCR+1) = IWRK
                       IWRK1 = IWRKT (ICRS1)
                       XWRK1 = XDATT (IWRK1)
                    Else
                       If (ICRS2 < NORD) Then
                          IWRKT (ICRS1) = IWRKT (ICRS)
                          IWRK1 = IWRKT (ICRS1)
                          XWRK1 = XDATT (IWRK1)
                       Endif
                    End If
                    ICRS1 = MIN (NORD, ICRS1 + 1)
                    ICRS2 = MIN (NORD, ICRS2 + 1)
                 End Do
              Else
                 Do ICRS = ISTRT (IDON1), IDON + 6
                    If (XDATT (IWRKT (ICRS)) >= XWRK1) Exit
                    IWRK = IWRKT (ICRS)
                    XWRK = XDATT (IWRK)
                    Do IDCR = ICRS1 - 1, 1, - 1
                       If (XWRK >= XDATT (IWRKT (IDCR))) Exit
                       IWRKT  (IDCR+1) = IWRKT (IDCR)
                    End Do
                    IWRKT (IDCR+1) = IWRK
                    IWRK1 = IWRKT (ICRS1)
                    XWRK1 = XDATT (IWRK1)
                 End Do
              End If
           End Do
           ires_med = IWRK1
           Return
        Else
           ires_med = IMED7
           Return
        End If
     Else                       !      If (NLEQ > NMED)
!                                          Not enough high values
        XWRK1 = -XHUGE
        NORD = NLEQ - NEQU - NMED + 1
        IDON1 = 0
        ICRS1 = 1
        ICRS2 = 0
        Do IDON = 1, NTRI, 7
           IDON1 = IDON1 + 1
           If (ICRS2 < NORD) Then
!
              Do ICRS = IDON, IENDT (IDON1)
                 If (XDATT(IWRKT (ICRS)) > XWRK1) Then
                    IWRK = IWRKT (ICRS)
                    XWRK = XDATT (IWRK)
                    IDCR = ICRS1 - 1
                    Do IDCR = ICRS1 - 1, 1, - 1
                       If (XWRK <= XDATT(IWRKT (IDCR))) Exit
                       IWRKT (IDCR+1) = IWRKT (IDCR)
                    End Do
                    IWRKT (IDCR+1) = IWRK
                    IWRK1 = IWRKT(ICRS1)
                    XWRK1 = XDATT(IWRK1)
                 Else
                    If (ICRS2 < NORD) Then
                       IWRKT (ICRS1) = IWRKT (ICRS)
                       IWRK1 = IWRKT(ICRS1)
                       XWRK1 = XDATT(IWRK1)
                    End If
                 End If
                 ICRS1 = MIN (NORD, ICRS1 + 1)
                 ICRS2 = MIN (NORD, ICRS2 + 1)
              End Do
           Else
              Do ICRS = IENDT (IDON1), IDON, -1
                 If (XDATT(IWRKT (ICRS)) <= XWRK1) Exit
                 IWRK = IWRKT (ICRS)
                 XWRK = XDATT (IWRK)
                 IDCR = ICRS1 - 1
                 Do IDCR = ICRS1 - 1, 1, - 1
                    If (XWRK <= XDATT(IWRKT (IDCR))) Exit
                    IWRKT (IDCR+1) = IWRKT (IDCR)
                 End Do
                 IWRKT (IDCR+1) = IWRK
                 IWRK1 = IWRKT(ICRS1)
                 XWRK1 = XDATT(IWRK1)
              End Do
           Endif
        End Do
!
        ires_med = IWRK1
        Return
     End If
