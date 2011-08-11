!                    ***************
                     FUNCTION QBBJ78
!                    ***************
!
     &( B     , IQBBJ )
!
!***********************************************************************
! TOMAWAC   V6P1                                   23/06/2011
!***********************************************************************
!
!brief    COMPUTES THE FRACTION OF BREAKING WAVES: QB.
!+                QB IS USED IN BATTJES AND JANSSEN (1978).
!
!reference  BATTJES AND JANSSEN (1978) :
!+                     "ENERGY LOSS AND SET-UP DUE TO BREAKING
!+                      OF RANDOM WAVES". ICCE'78.
!reference DINGEMANS (1983) :
!+                     "VERIFICATION OF NUMERICAL WAVE PROPAGATION
!+                      MODELS WITH FIELD MEASUREMENTS. CREDIZ
!+                      VERIFICATION HARINGVLIET".
!
!history  F.BECQ; M. BENOIT (EDF/DER/LNH)
!+        14/02/96
!+        V1P1
!+
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        13/07/2010
!+        V6P0
!+   Translation of French comments within the FORTRAN sources into
!+   English comments
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        21/08/2010
!+        V6P0
!+   Creation of DOXYGEN tags for automated documentation and
!+   cross-referencing of the FORTRAN sources
!
!history  G.MATTAROLO (EDF - LNHE)
!+        23/06/2011
!+        V6P1
!+   Translation of French names of the variables in argument
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| B              |-->| (H_RMS/H_MAX)
!| IQBBJ          |-->| INDEX OF THE SLECTED COMPUTATION METHOD
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      IMPLICIT NONE
!
!.....VARIABLES IN ARGUMENT
!     """"""""""""""""""""
      DOUBLE PRECISION QBBJ78, B
      INTEGER  IQBBJ
!
!.....LOCAL VARIABLES
!     """""""""""""""""
      DOUBLE PRECISION F     , CB    , EPS   , QMAX  , QMIN  , Q0
      DOUBLE PRECISION B2    , EXPO
!
!
      EPS = 1.D-7
!
      IF (B.GE.1.D0) THEN
        QBBJ78 = 1.D0
        RETURN
      ENDIF
!
      IF(IQBBJ.EQ.0) THEN
!       =========================
!       SOLVES BY DICHOTOMY
!       =========================
        QMIN  = 0.D0
        QMAX  = B
   10   CONTINUE
        QBBJ78 = (QMIN+QMAX)/2.D0
        F      = 1.D0 - QBBJ78 + B*B*DLOG(QBBJ78)
        IF (ABS(F).LT.EPS) RETURN
        IF (F.GT.0.D0) THEN
           QMAX = QBBJ78
        ELSE
           QMIN = QBBJ78
        ENDIF
        GOTO 10
!
      ELSEIF(IQBBJ.EQ.1) THEN
!       ======================================================
!       EXPLICIT FORMULATION 1 (INSPIRED FROM CREDIZ VERSION-1)
!       ======================================================
      CB = 0.5D0
        IF (B.GE.CB) THEN
          QBBJ78 = ((B-CB)/(1.D0-CB))**2
        ELSE
          QBBJ78 = 0.D0
        ENDIF
!
      ELSEIF(IQBBJ.EQ.2) THEN
!       ======================================================
!       EXPLICIT FORMULATION 2 (INSPIRED FROM CREDIZ VERSION-2)
!       ======================================================
        CB = 0.3D0
        IF (B.LT.CB) THEN
          QBBJ78 = 0.D0
        ELSEIF (B.LT.0.5D0) THEN
          B2     = B**2
          EXPO   = DEXP(-1.D0/B2)
          QBBJ78 = B2*EXPO/(B2-EXPO)
        ELSEIF (B.LT.0.9D0) THEN
          Q0     = (2.D0*B-1.D0)**2
          B2     = B**2
          EXPO   = DEXP((Q0-1.D0)/B2)
          QBBJ78 = Q0 - B2*(Q0-EXPO)/(B2-EXPO)
        ELSE
          QBBJ78 = (2.D0*B-1.D0)**2
        ENDIF
!
      ELSEIF(IQBBJ.EQ.3) THEN
!       ======================================================
!       EXPLICIT FORMULATION 3 (INSPIRED FROM CREDIZ VERSION-3)
!       ======================================================
        QBBJ78 = 2.4D0*B**7
!
      ENDIF
!
      RETURN
      END
