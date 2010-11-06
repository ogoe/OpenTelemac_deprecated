C             *******************************
              DOUBLE PRECISION FUNCTION ATANC
C             *******************************
C
     *(A,B)
C
C***********************************************************************
C BIEF VERSION 5.1          12/07/95    E. DAVID (LHF) 76 33 42 08
C
C***********************************************************************
C
C     FONCTIONS :  CALCUL D'UN ARC TANGENTE ENTRE 0 ET 2 PI.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C º      NOM       ºMODEº                   ROLE
C º________________º____º______________________________________________
C º   A            º -->º ANGLE DONT ON CHERCHE L'ARC TANGENTE.
C º   B            º -->º ANGLE PROCHE DE A
C º________________º____º______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES :
C
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION, INTENT(IN) :: A,B
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION PISUR2,D
C
      INTRINSIC ACOS,ATAN
C
C---------------------------------------------------------------
C
      PISUR2 = 0.5D0 * ACOS(-1.D0)
      D = ATAN(A)
C
10    CONTINUE
      IF (D-B.LT.-0.5D0*PISUR2) THEN
        D = D + PISUR2
        GO TO 10
      ENDIF
20    CONTINUE
      IF (D-B.GT.0.5D0*PISUR2) THEN
        D = D - PISUR2
        GO TO 20
      ENDIF
C
      ATANC = D
C
C---------------------------------------------------------------
C
      RETURN
      END
 
