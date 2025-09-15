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
module num_rec
  use healpix_types
  use misc_utils
  implicit none

  integer(i4b), parameter :: IZERO = 0_i4b, IONE = 1_i4b, ITWO = 2_i4b
  private


#ifndef NO64BITS
  interface isort
     module procedure isort, isort_8
  end interface
  interface iindexx
     module procedure iindexx, iindexx_8
  end interface
  interface sort
     module procedure sort, sort_8
  end interface
  interface indexx
     module procedure indexx, indexx_8
  end interface
#endif
  public :: dsvbksb, dsvdcmp, isort, sort, iindexx, indexx, othpl, pythag


contains
! f pythag
! s SVD
! s IPSORT
! s IPSORT_8
! S XPSORT
! s OTHPL
! s dsvbksb
! s dsvdcmp
! s isort
! s isort_8
! s sort
! S iindexx
! s iindexx_8
! S indexx
!
! The routines PYTHAG, SVD and IPSORT were taken from SLATEC:
! http://www.netlib.org/slatec/src/<routine>.f
! The SLATEC routines are in the public domain.
!
! Changes by MR:
!  - processed with f2f90 sed script
!  - replaced REAL with REAL(DP)
!  - replaced xx.0E0 with xx.0_DP
!  - rewrote some loops to avoid shared CONTINUE statements and loop end
!    statements which are not CONTINUE
!  - removed calls to XERMSG; the error code is checked in the calling routine
!
! Changes by EH:
!  psort has been edited from slatec's IPSORT
!   IX --> XX
!   ITEMP --> XTEMP
! Jan 2008: replaced integer by integer(I4B)
!           0, 1, 2 -> IZERO, IONE, ITWO
!           int(x)  -> int(x, kind=i4b)
! Mar 2017: PYTHAG: made public, replaced xx.0 with xx.0_dp
!
! The routine OTHPL was taken from:
!   S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
!    http://jin.ece.uiuc.edu/routines/routines.html
! editions by EH:
! replaced implicit declaration by explicit ones,
! added intent of dummy arguments
! added tests (n>0)
! -----------------------------------------------------------
!
!DECK PYTHAG
      FUNCTION PYTHAG (A, B) result(pyth)
!***BEGIN PROLOGUE  PYTHAG
!***SUBSIDIARY
!***PURPOSE  Compute the complex square root of a complex number without
!            destructive overflow or underflow.
!***LIBRARY   SLATEC
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     Finds sqrt(A**2+B**2) without overflow or destructive underflow
!
!***SEE ALSO  EISDOC
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   811101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  PYTHAG
      REAL(DP):: PYTH
      REAL(DP), intent(IN):: A,B
!
      REAL(DP):: P,Q,R,S,T
!***FIRST EXECUTABLE STATEMENT  PYTHAG
      P = MAX(ABS(A),ABS(B))
      Q = MIN(ABS(A),ABS(B))
      IF (Q .EQ. 0.0_dp) GO TO 20
   10 CONTINUE
         R = (Q/P)**ITWO
         T = 4.0_dp + R
         IF (T .EQ. 4.0_dp) GO TO 20
         S = R/T
         P = P + 2.0_dp*P*S
         Q = Q*S
      GO TO 10
   20 PYTH = P
      RETURN
      END FUNCTION PYTHAG

!DECK SVD
      SUBROUTINE SVD (NM, M, N, A, W, MATU, U, MATV, V, IERR, RV1)
!***BEGIN PROLOGUE  SVD
!***SUBSIDIARY
!***PURPOSE  Perform the singular value decomposition of a rectangular
!            matrix.
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (SVD-S)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure SVD,
!     NUM. MATH. 14, 403-420(1970) by Golub and Reinsch.
!     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 134-151(1971).
!
!     This subroutine determines the singular value decomposition
!          T
!     A=USV  of a REAL M by N rectangular matrix.  Householder
!     bidiagonalization and a variant of the QR algorithm are used.
!
!     On Input
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, A, U and V, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!          Note that NM must be at least as large as the maximum
!          of M and N.
!
!        M is the number of rows of A and U.
!
!        N is the number of columns of A and U and the order of V.
!
!        A contains the rectangular input matrix to be decomposed.  A is
!          a two-dimensional REAL array, dimensioned A(NM,N).
!
!        MATU should be set to .TRUE. if the U matrix in the
!          decomposition is desired, and to .FALSE. otherwise.
!          MATU is a LOGICAL variable.
!
!        MATV should be set to .TRUE. if the V matrix in the
!          decomposition is desired, and to .FALSE. otherwise.
!          MATV is a LOGICAL variable.
!
!     On Output
!
!        A is unaltered (unless overwritten by U or V).
!
!        W contains the N (non-negative) singular values of A (the
!          diagonal elements of S).  They are unordered.  If an
!          error exit is made, the singular values should be correct
!          for indices IERR+1, IERR+2, ..., N.  W is a one-dimensional
!          REAL array, dimensioned W(N).
!
!        U contains the matrix U (orthogonal column vectors) of the
!          decomposition if MATU has been set to .TRUE.  Otherwise,
!          U is used as a temporary array.  U may coincide with A.
!          If an error exit is made, the columns of U corresponding
!          to indices of correct singular values should be correct.
!          U is a two-dimensional REAL array, dimensioned U(NM,N).
!
!        V contains the matrix V (orthogonal) of the decomposition if
!          MATV has been set to .TRUE.  Otherwise, V is not referenced.
!          V may also coincide with A if U does not.  If an error
!          exit is made, the columns of V corresponding to indices of
!          correct singular values should be correct.  V is a two-
!          dimensional REAL array, dimensioned V(NM,N).
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          K          if the K-th singular value has not been
!                     determined after 30 iterations.
!
!        RV1 is a one-dimensional REAL array used for temporary storage,
!          dimensioned RV1(N).
!
!     CALLS PYTHAG(A,B) for sqrt(A**2 + B**2).
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***SEE ALSO  EISDOC
!***ROUTINES CALLED  PYTHAG
!***REVISION HISTORY  (YYMMDD)
!   811101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  SVD
!
      INTEGER(I4B) :: I,J,K,L,M,N,II,I1,KK,K1,LL,L1,MN,NM,ITS,IERR
      REAL(DP) A(NM,*),W(*),U(NM,*),V(NM,*),RV1(*)
      REAL(DP) C,F,G,H,S,X,Y,Z,SCALE,S1
!      REAL(DP) PYTHAG
      LOGICAL MATU,MATV
!
!***FIRST EXECUTABLE STATEMENT  SVD
      IERR = 0
!
      DO I = IONE, M
!
         DO J = IONE, N
            U(I,J) = A(I,J)
         END DO
      END DO
!     .......... HOUSEHOLDER REDUCTION TO BIDIAGONAL FORM ..........
      G = 0.0_DP
      SCALE = 0.0_DP
      S1 = 0.0_DP
!
      DO 300 I = IONE, N
         L = I + IONE
         RV1(I) = SCALE * G
         G = 0.0_DP
         S = 0.0_DP
         SCALE = 0.0_DP
         IF (I .GT. M) GO TO 210
!
         DO K = I, M
           SCALE = SCALE + ABS(U(K,I))
         END DO
!
         IF (SCALE .EQ. 0.0_DP) GO TO 210
!
         DO 130 K = I, M
            U(K,I) = U(K,I) / SCALE
            S = S + U(K,I)**2
  130    CONTINUE
!
         F = U(I,I)
         G = -SIGN(SQRT(S),F)
         H = F * G - S
         U(I,I) = F - G
         IF (I .EQ. N) GO TO 190
!
         DO 150 J = L, N
            S = 0.0_DP
!
            DO K = I, M
              S = S + U(K,I) * U(K,J)
            END DO
!
            F = S / H
!
            DO K = I, M
              U(K,J) = U(K,J) + F * U(K,I)
            END DO
  150    CONTINUE
!
  190    DO K = I, M
           U(K,I) = SCALE * U(K,I)
         END DO
!
  210    W(I) = SCALE * G
         G = 0.0_DP
         S = 0.0_DP
         SCALE = 0.0_DP
         IF (I .GT. M .OR. I .EQ. N) GO TO 290
!
         DO K = L, N
           SCALE = SCALE + ABS(U(I,K))
         END DO
!
         IF (SCALE .EQ. 0.0_DP) GO TO 290
!
         DO 230 K = L, N
            U(I,K) = U(I,K) / SCALE
            S = S + U(I,K)**2
  230    CONTINUE
!
         F = U(I,L)
         G = -SIGN(SQRT(S),F)
         H = F * G - S
         U(I,L) = F - G
!
         DO K = L, N
           RV1(K) = U(I,K) / H
         END DO
!
         IF (I .EQ. M) GO TO 270
!
         DO 260 J = L, M
            S = 0.0_DP
!
            DO K = L, N
              S = S + U(J,K) * U(I,K)
            END DO
!
            DO K = L, N
               U(J,K) = U(J,K) + S * RV1(K)
            END DO
  260    CONTINUE
!
  270    DO K = L, N
           U(I,K) = SCALE * U(I,K)
         END DO
!
  290    S1 = MAX(S1,ABS(W(I))+ABS(RV1(I)))
  300 CONTINUE
!     .......... ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS ..........
      IF (.NOT. MATV) GO TO 410
!     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 400 II = IONE, N
         I = N + IONE - II
         IF (I .EQ. N) GO TO 390
         IF (G .EQ. 0.0_DP) GO TO 360
!
         DO J = L, N
!     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
           V(J,I) = (U(I,J) / U(I,L)) / G
         END DO
!
         DO 350 J = L, N
            S = 0.0_DP
!
            DO K = L, N
              S = S + U(I,K) * V(K,J)
            END DO
!
            DO K = L, N
               V(K,J) = V(K,J) + S * V(K,I)
            END DO
  350    CONTINUE
!
  360    DO 380 J = L, N
            V(I,J) = 0.0_DP
            V(J,I) = 0.0_DP
  380    CONTINUE
!
  390    V(I,I) = 1.0_DP
         G = RV1(I)
         L = I
  400 CONTINUE
!     .......... ACCUMULATION OF LEFT-HAND TRANSFORMATIONS ..........
  410 IF (.NOT. MATU) GO TO 510
!     ..........FOR I=MIN(M,N) STEP -1 UNTIL 1 DO -- ..........
      MN = N
      IF (M .LT. N) MN = M
!
      DO 500 II = IONE, MN
         I = MN + IONE - II
         L = I + IONE
         G = W(I)
         IF (I .EQ. N) GO TO 430
!
         DO J = L, N
           U(I,J) = 0.0_DP
         END DO
!
  430    IF (G .EQ. 0.0_DP) GO TO 475
         IF (I .EQ. MN) GO TO 460
!
         DO 450 J = L, N
            S = 0.0_DP
!
            DO K = L, M
              S = S + U(K,I) * U(K,J)
            END DO
!     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
            F = (S / U(I,I)) / G
!
            DO K = I, M
               U(K,J) = U(K,J) + F * U(K,I)
            END DO
  450    CONTINUE
!
  460    DO J = I, M
           U(J,I) = U(J,I) / G
         END DO
!
         GO TO 490
!
  475    DO J = I, M
           U(J,I) = 0.0_DP
         END DO
!
  490    U(I,I) = U(I,I) + 1.0_DP
  500 CONTINUE
!     .......... DIAGONALIZATION OF THE BIDIAGONAL FORM ..........
  510 CONTINUE
!     .......... FOR K=N STEP -1 UNTIL 1 DO -- ..........
      DO 700 KK = IONE, N
         K1 = N - KK
         K = K1 + IONE
         ITS = 0
!     .......... TEST FOR SPLITTING.
!                FOR L=K STEP -1 UNTIL 1 DO -- ..........
  520    DO 530 LL = IONE, K
            L1 = K - LL
            L = L1 + IONE
            IF (S1 + ABS(RV1(L)) .EQ. S1) GO TO 565
!     .......... RV1(1) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP ..........
            IF (S1 + ABS(W(L1)) .EQ. S1) GO TO 540
  530    CONTINUE
!     .......... CANCELLATION OF RV1(L) IF L GREATER THAN 1 ..........
  540    C = 0.0_DP
         S = 1.0_DP
!
         DO 560 I = L, K
            F = S * RV1(I)
            RV1(I) = C * RV1(I)
            IF (S1 + ABS(F) .EQ. S1) GO TO 565
            G = W(I)
            H = PYTHAG(F,G)
            W(I) = H
            C = G / H
            S = -F / H
            IF (.NOT. MATU) GO TO 560
!
            DO 550 J = IONE, M
               Y = U(J,L1)
               Z = U(J,I)
               U(J,L1) = Y * C + Z * S
               U(J,I) = -Y * S + Z * C
  550       CONTINUE
!
  560    CONTINUE
!     .......... TEST FOR CONVERGENCE ..........
  565    Z = W(K)
         IF (L .EQ. K) GO TO 650
!     .......... SHIFT FROM BOTTOM 2 BY 2 MINOR ..........
         IF (ITS .EQ. 30) GO TO 1000
         ITS = ITS + IONE
         X = W(L)
         Y = W(K1)
         G = RV1(K1)
         H = RV1(K)
         F = 0.5_DP * (((G + Z) / H) * ((G - Z) / Y) + Y / H - H / Y)
         G = PYTHAG(F,1.0_DP)
         F = X - (Z / X) * Z + (H / X) * (Y / (F + SIGN(G,F)) - H)
!     .......... NEXT QR TRANSFORMATION ..........
         C = 1.0_DP
         S = 1.0_DP
!
         DO 600 I1 = L, K1
            I = I1 + IONE
            G = RV1(I)
            Y = W(I)
            H = S * G
            G = C * G
            Z = PYTHAG(F,H)
            RV1(I1) = Z
            C = F / Z
            S = H / Z
            F = X * C + G * S
            G = -X * S + G * C
            H = Y * S
            Y = Y * C
            IF (.NOT. MATV) GO TO 575
!
            DO 570 J = IONE, N
               X = V(J,I1)
               Z = V(J,I)
               V(J,I1) = X * C + Z * S
               V(J,I) = -X * S + Z * C
  570       CONTINUE
!
  575       Z = PYTHAG(F,H)
            W(I1) = Z
!     .......... ROTATION CAN BE ARBITRARY IF Z IS ZERO ..........
            IF (Z .EQ. 0.0_DP) GO TO 580
            C = F / Z
            S = H / Z
  580       F = C * G + S * Y
            X = -S * G + C * Y
            IF (.NOT. MATU) GO TO 600
!
            DO 590 J = IONE, M
               Y = U(J,I1)
               Z = U(J,I)
               U(J,I1) = Y * C + Z * S
               U(J,I) = -Y * S + Z * C
  590       CONTINUE
!
  600    CONTINUE
!
         RV1(L) = 0.0_DP
         RV1(K) = F
         W(K) = X
         GO TO 520
!     .......... CONVERGENCE ..........
  650    IF (Z .GE. 0.0_DP) GO TO 700
!     .......... W(K) IS MADE NON-NEGATIVE ..........
         W(K) = -Z
         IF (.NOT. MATV) GO TO 700
!
         DO J = IONE, N
           V(J,K) = -V(J,K)
         END DO
!
  700 CONTINUE
!
      GO TO 1001
!     .......... SET ERROR -- NO CONVERGENCE TO A
!                SINGULAR VALUE AFTER 30 ITERATIONS ..........
 1000 IERR = K
 1001 RETURN
      END SUBROUTINE SVD


!=====================================
!   ipsort   (ix, n, iperm, kflag, ier)
!   ipsort_8 (ix, n, iperm, kflag, ier)
!   xpsort   (ix, n, iperm, kflag, ier)
!   xpsort_8 (ix, n, iperm, kflag, ier)
!
!***PURPOSE  Return the permutation vector generated by sorting a given
!            array and, optionally, rearrange the elements of the array.
!            The array may be sorted in increasing or decreasing order.
!            A slightly modified quicksort algorithm is used.
!***LIBRARY   SLATEC
!***CATEGORY  N6A1A, N6A2A
!***TYPE      INTEGER (SPSORT-S, DPSORT-D, IPSORT-I, HPSORT-H)
!***KEYWORDS  NUMBER SORTING, PASSIVE SORTING, SINGLETON QUICKSORT, SORT
!***AUTHOR  Jones, R. E., (SNLA)
!           Kahaner, D. K., (NBS)
!           Rhoads, G. S., (NBS)
!           Wisniewski, J. A., (SNLA)
!***DESCRIPTION
!
!   IPSORT returns the permutation vector IPERM generated by sorting
!   the array IX and, optionally, rearranges the values in IX.  IX may
!   be sorted in increasing or decreasing order.  A slightly modified
!   quicksort algorithm is used.
!
!   IPERM is such that IX(IPERM(I)) is the Ith value in the
!   rearrangement of IX.  IPERM may be applied to another array by
!   calling IPPERM, SPPERM, DPPERM or HPPERM.
!
!   The main difference between IPSORT and its active sorting equivalent
!   ISORT is that the data are referenced indirectly rather than
!   directly.  Therefore, IPSORT should require approximately twice as
!   long to execute as ISORT.  However, IPSORT is more general.
!
!   Description of Parameters
!      IX - input/output -- integer array of values to be sorted.
!           If ABS(KFLAG) = 2, then the values in IX will be
!           rearranged on output; otherwise, they are unchanged.
!      N  - input -- number of values in array IX to be sorted.
!      IPERM - output -- permutation array such that IPERM(I) is the
!              index of the value in the original order of the
!              IX array that is in the Ith location in the sorted
!              order.
!      KFLAG - input -- control parameter:
!            =  2  means return the permutation vector resulting from
!                  sorting IX in increasing order and sort IX also.
!            =  1  means return the permutation vector resulting from
!                  sorting IX in increasing order and do not sort IX.
!            = -1  means return the permutation vector resulting from
!                  sorting IX in decreasing order and do not sort IX.
!            = -2  means return the permutation vector resulting from
!                  sorting IX in decreasing order and sort IX also.
!      IER - output -- error indicator:
!          =  0  if no error,
!          =  1  if N is zero or negative,
!          =  2  if KFLAG is not 2, 1, -1, or -2.
!***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
!                 for sorting with minimal storage, Communications of
!                 the ACM, 12, 3 (1969), pp. 185-187.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   761101  DATE WRITTEN
!   761118  Modified by John A. Wisniewski to use the Singleton
!           quicksort algorithm.
!   810801  Further modified by David K. Kahaner.
!   870423  Modified by Gregory S. Rhoads for passive sorting with the
!           option for the rearrangement of the original data.
!   890620  Algorithm for rearranging the data vector corrected by R.
!           Boisvert.
!   890622  Prologue upgraded to Version 4.0 style by D. Lozier.
!   891128  Error when KFLAG.LT.0 and N=1 corrected by R. Boisvert.
!   920507  Modified by M. McClain to revise prologue text.
!   920818  Declarations section rebuilt and code restructured to use
!           IF-THEN-ELSE-ENDIF.  (SMR, WRB)
!
! 2009-06-17: uses generic source file :  source_gen_sort_inc.f90
!  the array to sort is now generically gen_xx
!=====================================
      SUBROUTINE IPSORT (gen_xx, N, IPERM, KFLAG, IER)
        ! integer 4 array with less than 2^31 elements to sort
!=====================================
      integer(i4b), parameter :: MKD = I4B
      INTEGER(MKD), intent(in)  :: N
      INTEGER(I4B), intent(in)  :: KFLAG
      INTEGER(I4B), intent(out) :: IER
      INTEGER(MKD),                dimension(1:*) :: gen_xx
      INTEGER(MKD), intent(out),   dimension(1:*) :: IPERM
!     .. Local Scalars ..
      REAL(DP) ::  R
      INTEGER(MKD) :: gen_temp
      INTEGER(MKD) :: I, IJ, INDX, INDX0, ISTRT, J, K, KK, L, LM, LMT, M, NN
!     .. Local Arrays ..
      INTEGER(MKD) :: IL(21), IU(21)

#include "source_gen_sort_inc.f90"

      RETURN
      END SUBROUTINE IPSORT

!=====================================
      SUBROUTINE IPSORT_8 (gen_xx, N, IPERM, KFLAG, IER)
        ! integer 8 array with any number of elements to sort
!=====================================
      integer(i4b), parameter :: MKD = I8B
      INTEGER(MKD), intent(in)  :: N
      INTEGER(I4B), intent(in)  :: KFLAG
      INTEGER(I4B), intent(out) :: IER
      INTEGER(MKD),                dimension(1:*) :: gen_xx
      INTEGER(MKD), intent(out),   dimension(1:*) :: IPERM
!     .. Local Scalars ..
      REAL(DP) ::  R
      INTEGER(MKD) :: gen_temp
      INTEGER(MKD) :: I, IJ, INDX, INDX0, ISTRT, J, K, KK, L, LM, LMT, M, NN
!     .. Local Arrays ..
      INTEGER(MKD) :: IL(21), IU(21)

#include "source_gen_sort_inc.f90"

      RETURN
      END SUBROUTINE IPSORT_8

!=====================================
      SUBROUTINE XPSORT (gen_xx, N, IPERM, KFLAG, IER)
        ! double precision real with less than 2^31 elements to sort
! ========================================================
        integer(i4b), parameter :: MKD = I4B
      INTEGER(MKD), intent(in)  :: N
      INTEGER(I4B), intent(in)  :: KFLAG
      INTEGER(I4B), intent(out) :: IER
      REAL(DP),                    dimension(1:*) :: gen_xx
      INTEGER(MKD), intent(out),   dimension(1:*) :: IPERM
!     .. Local Scalars ..
      REAL(DP) :: R
      real(DP) :: gen_temp
      INTEGER(I4B) :: I, IJ, INDX, INDX0, ISTRT, J, K, KK, L, LM, LMT, M, NN
!     .. Local Arrays ..
      INTEGER(I4B) :: IL(21), IU(21)

#include "source_gen_sort_inc.f90"

      RETURN
    END SUBROUTINE XPSORT
!=========================================

!=====================================
      SUBROUTINE XPSORT_8 (gen_xx, N, IPERM, KFLAG, IER)
        ! double precision real with any number of elements to sort
! ========================================================
      integer(i4b), parameter :: MKD = I8B
      INTEGER(MKD), intent(in)  :: N
      INTEGER(I4B), intent(in)  :: KFLAG
      INTEGER(I4B), intent(out) :: IER
      REAL(DP),                    dimension(1:*) :: gen_xx
      INTEGER(MKD), intent(out),   dimension(1:*) :: IPERM
!     .. Local Scalars ..
      REAL(DP) :: R
      real(DP) :: gen_temp
      INTEGER(I4B) :: I, IJ, INDX, INDX0, ISTRT, J, K, KK, L, LM, LMT, M, NN
!     .. Local Arrays ..
      INTEGER(I4B) :: IL(21), IU(21)

#include "source_gen_sort_inc.f90"

      RETURN
    END SUBROUTINE XPSORT_8
!=========================================


    SUBROUTINE OTHPL(KF,N,X,PL,DPL)
      !
      !       ==========================================================
      !       Purpose: Compute orthogonal polynomials: Tn(x) or Un(x),
      !                or Ln(x) or Hn(x), and their derivatives
      !       Input :  KF --- Function code
      !                       KF=1 for Chebyshev polynomial Tn(x)
      !                       KF=2 for Chebyshev polynomial Un(x)
      !                       KF=3 for Laguerre polynomial Ln(x)
      !                       KF=4 for Hermite polynomial Hn(x)
      !                n ---  Order of orthogonal polynomials
      !                x ---  Argument of orthogonal polynomials
      !       Output:  PL(n) --- Tn(x) or Un(x) or Ln(x) or Hn(x)
      !                DPL(n)--- Tn'(x) or Un'(x) or Ln'(x) or Hn'(x)
      !       =========================================================
      ! copyright: 
      !   S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
      !    http://jin.ece.uiuc.edu/routines/routines.html
      !--------------------------------------------------------------------------
      use healpix_types
      IMPLICIT none
      integer(i4b), intent(in) :: kf, n
      real(DP), intent(in) :: x
      real(DP), dimension(0:), intent(out) :: pl(0:),dpl(0:)
      !
      real(DP) :: a, b, c, y0, y1, dy0, dy1, yn, dyn
      integer(i4b) :: k

      a = 2.0d0
      b = 0.0d0
      c = 1.0d0
      y0 = 1.0d0
      y1 = 2.0d0*x
      dy0 = 0.0d0
      dy1 = 2.0d0
      pl(0) = 1.0d0
      dpl(0) = 0.0d0
      if (n > 0) then
         pl(1) = 2.0d0*x
         dpl(1) = 2.0d0
      endif
      if (kf == 1) then
         y1 = x
         dy1 = 1.0d0
         if (n > 0) then
            pl(1) = x
            dpl(1) = 1.0d0
         endif
      else if (kf == 3) then
         y1 = 1.0d0-x
         dy1 = -1.0d0
         if (n > 0) then
            pl(1) = 1.0d0-x
            dpl(1) = -1.0d0
         endif
      endif
      do k = 2,n
         if (kf == 3) then
            a = -1.0d0/k
            b = 2.0d0+a
            c = 1.0d0+a
         else if (kf == 4) then
            c = 2.0d0*(k-1.0d0)
         endif
         yn = (a*x+b)*y1-c*y0
         dyn = a*y1+(a*x+b)*dy1-c*dy0
         pl(k) = yn
         dpl(k) = dyn
         y0 = y1
         y1 = yn
         dy0 = dy1
         dy1 = dyn
      enddo
      return

    END SUBROUTINE OTHPL

!=========================================

  subroutine dsvbksb(u,w,v,m,n,mp,np,b,x)
    integer(i4b), intent(in) :: m,mp,n,np
    real(dp), intent(in) ::  b(mp),u(mp,np),v(np,np),w(np)
    real(dp), intent(out) :: x(np)
    integer(i4b) ::  i,j,jj
    real(dp) ::  s,tmp(n)
    do j=IONE,n
       s = 0d0
       if (w(j)/=0d0) then
          do i=IONE,m
             s=s+u(i,j)*b(i)
          end do
          s=s/w(j)
       endif
       tmp(j)=s
    end do
    do j=IONE,n
       s=0.0d0
       do jj=IONE,n
          s=s+v(j,jj)*tmp(jj)
       end do
       x(j)=s
    end do
  end subroutine dsvbksb

  subroutine dsvdcmp(a,m,n,mp,np,w,v)
    integer(i4b), intent(in):: m,mp,n,np
    real(dp), intent(inout) ::  a(mp,np),v(np,np),w(np)
    integer(i4b) :: ierr
    real(dp) rv1(np),u(mp,np)

    call svd(mp,m,n,a,w,.true.,u,.true.,v,ierr,rv1)
    call assert(ierr==0,"error in svd()")
    a=u
  end subroutine dsvdcmp

  subroutine isort(n,ra)
    integer(i4b), intent(in) :: n
    integer(i4b), intent(inout), dimension(1:) :: ra
    integer(i4b) :: indx(n)
    integer(i4b) :: ier
    call ipsort (ra,n,indx,ITWO,ier)
    call assert (ier==IZERO, "error in ipsort()")
  end subroutine isort

  subroutine isort_8(n,ra)
    integer(i8b), intent(in) :: n
    integer(i8b), intent(inout), dimension(1:) :: ra
    integer(i8b) :: indx(n)
    integer(i4b) :: ier
    call ipsort_8(ra,n,indx,ITWO,ier)
    call assert (ier==IZERO, "error in ipsort_8()")
  end subroutine isort_8

  subroutine sort(n,ra)
    integer(i4b), intent(in) :: n
    real(dp), intent(inout), dimension(1:) :: ra
    integer(i4b) :: indx(n)
    integer(i4b) :: ier
    call xpsort (ra,n,indx,ITWO,ier)
    call assert (ier==IZERO, "error in xpsort()")
  end subroutine sort

  subroutine sort_8(n,ra)
    integer(i8b), intent(in) :: n
    real(dp), intent(inout), dimension(1:) :: ra
    integer(i8b) :: indx(n)
    integer(i4b) :: ier
    call xpsort_8(ra,n,indx,ITWO,ier)
    call assert (ier==IZERO, "error in xpsort_8()")
  end subroutine sort_8

  SUBROUTINE iindexx(n,arr,indx)
    INTEGER(I4B),                 intent(IN)  :: n
    INTEGER(I4B), dimension(1:n), intent(IN)  :: arr
    integer(i4b), dimension(1:n), intent(out) :: indx
    integer(i4b) :: ier
    call ipsort (arr,n,indx,IONE,ier)
    call assert (ier==IZERO, "error in ipsort()")
  END SUBROUTINE iindexx

  SUBROUTINE iindexx_8(n,arr,indx)
    INTEGER(I8B),                 intent(IN)  :: n
    INTEGER(I8B), dimension(1:n), intent(IN)  :: arr
    integer(i8b), dimension(1:n), intent(out) :: indx
    integer(i4b) :: ier
    call ipsort_8(arr,n,indx,IONE,ier)
    call assert (ier==IZERO, "error in ipsort_8()")
  END SUBROUTINE iindexx_8

  SUBROUTINE indexx(n,arr,indx)
    INTEGER(I4B),                  intent(IN)  :: n
    real(DP),      dimension(1:n), intent(IN)  :: arr
    integer(i4b),  dimension(1:n), intent(out) :: indx
    integer(i4b) ::  ier
    call xpsort (arr,n,indx,IONE,ier)
    call assert (ier==IZERO, "error in xpsort()")
  END SUBROUTINE indexx

  SUBROUTINE indexx_8(n,arr,indx)
    INTEGER(I8B),                  intent(IN)  :: n
    real(DP),      dimension(1:n), intent(IN)  :: arr
    integer(i8b),  dimension(1:n), intent(out) :: indx
    integer(i4b) ::  ier
    call xpsort_8(arr,n,indx,IONE,ier)
    call assert  (ier==IZERO, "error in xpsort()")
  END SUBROUTINE indexx_8

end module num_rec
