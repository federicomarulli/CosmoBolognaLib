!-----------------------------------------------------------------------------
!
!  Copyright (C) 1997-2013 Krzysztof M. Gorski, Eric Hivon,
!                          Benjamin D. Wandelt, Anthony J. Banday, 
!                          Matthias Bartelmann, Hans K. Eriksen, 
!                          Frode K. Hansen, Martin Reinecke
!
!
!  This file is part of HEALPix.
!
!  HEALPix is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  HEALPix is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with HEALPix; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!  For more information about HEALPix see http://healpix.sourceforge.net
!
!-----------------------------------------------------------------------------
!
! generic source for ipsort, xpsort, ipsort_8, ....
! the array to sort is gen_xx and can be integer, long integer, float, double ...
!
      IER = IZERO
      NN = N
      IF (NN .LT. IONE) THEN
         IER = IONE
         RETURN
      ENDIF
      KK = ABS(KFLAG)
      IF (KK.NE.IONE .AND. KK.NE.2) THEN
         IER = ITWO
         RETURN
      ENDIF
!
!     Initialize permutation vector
!
      DO 10 I=IONE,NN
         IPERM(I) = I
   10 CONTINUE
!
!     Return if only one value is to be sorted
!
      IF (NN .EQ. IONE) RETURN
!
!     Alter array GEN_XX to get decreasing order if needed
!
      IF (KFLAG .LE. -IONE) THEN
         DO 20 I=IONE,NN
            GEN_XX(I) = -GEN_XX(I)
   20    CONTINUE
      ENDIF
!
!     Sort GEN_XX only
!
      M = IONE
      I = IONE
      J = NN
      R = .375_DP  ! 15/40
!
   30 IF (I .EQ. J) GO TO 80
      IF (R .LE. 0.5898437_DP) THEN
         R = R+3.90625E-2  ! R = R + 5/128
      ELSE
         R = R-0.21875_DP  ! R = R - 7/32
      ENDIF
!
   40 K = I
!
!     Select a central element of the array and save it in location L
!
      IJ = I + INT((J-I)*R, kind=MKD)
      LM = IPERM(IJ)
!
!     If first element of array is greater than LM, interchange with LM
!
      IF (GEN_XX(IPERM(I)) .GT. GEN_XX(LM)) THEN
         IPERM(IJ) = IPERM(I)
         IPERM(I) = LM
         LM = IPERM(IJ)
      ENDIF
      L = J
!
!     If last element of array is less than LM, interchange with LM
!
      IF (GEN_XX(IPERM(J)) .LT. GEN_XX(LM)) THEN
         IPERM(IJ) = IPERM(J)
         IPERM(J) = LM
         LM = IPERM(IJ)
!
!        If first element of array is greater than LM, interchange
!        with LM
!
         IF (GEN_XX(IPERM(I)) .GT. GEN_XX(LM)) THEN
            IPERM(IJ) = IPERM(I)
            IPERM(I) = LM
            LM = IPERM(IJ)
         ENDIF
      ENDIF
      GO TO 60
   50 LMT = IPERM(L)
      IPERM(L) = IPERM(K)
      IPERM(K) = LMT
!
!     Find an element in the second half of the array which is smaller
!     than LM
!
   60 L = L-IONE
      IF (GEN_XX(IPERM(L)) .GT. GEN_XX(LM)) GO TO 60
!
!     Find an element in the first half of the array which is greater
!     than LM
!
   70 K = K+IONE
      IF (GEN_XX(IPERM(K)) .LT. GEN_XX(LM)) GO TO 70
!
!     Interchange these elements
!
      IF (K .LE. L) GO TO 50
!
!     Save upper and lower subscripts of the array yet to be sorted
!
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+IONE
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+IONE
      ENDIF
      GO TO 90
!
!     Begin again on another portion of the unsorted array
!
   80 M = M-IONE
      IF (M .EQ. IZERO) GO TO 120
      I = IL(M)
      J = IU(M)
!
   90 IF (J-I .GE. IONE) GO TO 40
      IF (I .EQ. IONE) GO TO 30
      I = I-IONE
!
  100 I = I+IONE
      IF (I .EQ. J) GO TO 80
      LM = IPERM(I+IONE)
      IF (GEN_XX(IPERM(I)) .LE. GEN_XX(LM)) GO TO 100
      K = I
!
  110 IPERM(K+IONE) = IPERM(K)
      K = K-IONE
!
      IF (GEN_XX(LM) .LT. GEN_XX(IPERM(K))) GO TO 110
      IPERM(K+IONE) = LM
      GO TO 100
!
!     Clean up
!
  120 IF (KFLAG .LE. -IONE) THEN
         DO 130 I=IONE,NN
            GEN_XX(I) = -GEN_XX(I)
  130    CONTINUE
      ENDIF
!
!     Rearrange the values of GEN_XX if desired
!
      IF (KK .EQ. ITWO) THEN
!
!        Use the IPERM vector as a flag.
!        If IPERM(I) < 0, then the I-th value is in correct location
!
         DO 150 ISTRT=IONE,NN
            IF (IPERM(ISTRT) .GE. IZERO)  THEN
               INDX = ISTRT
               INDX0 = INDX
               GEN_TEMP = GEN_XX(ISTRT)
  140          IF (IPERM(INDX) .GT. IZERO)  THEN
                  GEN_XX(INDX) = GEN_XX(IPERM(INDX))
                  INDX0 = INDX
                  IPERM(INDX) = -IPERM(INDX)
                  INDX = ABS(IPERM(INDX))
                  GO TO 140
               ENDIF
               GEN_XX(INDX0) = GEN_TEMP
            ENDIF
  150    CONTINUE
!
!        Revert the signs of the IPERM values
!
         DO 160 I=IONE,NN
            IPERM(I) = -IPERM(I)
  160    CONTINUE
!
      ENDIF
!
