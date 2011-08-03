!                    *********************
                     INTEGER FUNCTION NEXT
!                    *********************
!
     &( ICOL , LIGNE )
!
!***********************************************************************
! DAMOCLES   V6P0                                   21/08/2010
!***********************************************************************
!
!brief    RETURNS THE INDEX OF THE 1ST NON-WHITE, NON-TABULATION
!+             AND NON-COMMENT CHARACTER OF THE LINE, STARTING FROM
!+             COLUMN ICOL.
!+             IF THERE ARE NONE, SCANS THE NEXT LINE
!+             IF CANNOT FIND ANY, NEXT = LONGLI + 1
!
!note     PORTABILITY : IBM,CRAY,HP,SUN
!
!history  O. QUIQUEMPOIX (LNH)
!+        15/12/1993
!+
!+
!
!history  J.M. HERVOUET (LNH); A. YESSAYAN
!+        16/08/1994
!+        V5P1
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
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| ICOL           |<->| POSITION COURANTE DU POINTEUR DANS LA LIGNE
!| LIGNE          |<->| LIGNE EN COURS DE DECODAGE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      IMPLICIT NONE
!
      INTEGER       ICOL
      CHARACTER*(*) LIGNE*(*)
!
      INTEGER       LNG,LU
      INTEGER       NLIGN,LONGLI
      INTEGER       NFIC
      LOGICAL       ERREUR,RETOUR
!
!-----------------------------------------------------------------------
!
      INTEGER       I,J
      CHARACTER*1   TABUL
!
!-----------------------------------------------------------------------
!
      COMMON / DCINFO / LNG,LU
      COMMON / DCRARE / ERREUR,RETOUR
      COMMON / DCMLIG / NLIGN,LONGLI
      COMMON / DCCHIE / NFIC
!
      INTRINSIC CHAR
!
!***********************************************************************
!                                    RCS AND SCCS MARKING
!
!***********************************************************************
!
      TABUL = CHAR(9)
      NEXT  = LONGLI + 1
      I     = ICOL -1
!
100   CONTINUE
      I = I + 1
!
      IF(I.GT.LONGLI) THEN
99        NLIGN = NLIGN + 1
          READ(NFIC,END=900,ERR=998,FMT='(A)') LIGNE
!         DOES NOT CONSIDER A LINE STARTING WITH '/'
          IF(LIGNE(1:1).EQ.'/') GO TO 99
          I = 0
          GO TO 100
      ENDIF
!
!
! CASE OF A WHITE OR TABULATION CHARACTER
      IF (LIGNE(I:I).EQ.' '.OR.LIGNE(I:I).EQ.TABUL) GOTO 100
!
!-----------------------------------------------------------------------
!
!          DOES NOT CONSIDER THE COMMENTED LINES:
!
!          IF ( LIGNE(I:I).NE.'/'.OR.I.GE.LONGLI ) THEN
           IF (LIGNE(I:I).NE.'/') THEN
                NEXT = I
                GO TO 1000
           ELSE
                DO 110 J = I+1 , LONGLI
                     IF ( LIGNE(J:J).EQ.'/' ) THEN
                          I = J
                          GO TO 100
                     ENDIF
110             CONTINUE
                I = LONGLI
                GO TO 100
            ENDIF
!
!-----------------------------------------------------------------------
!
998   CONTINUE
      IF(LNG.EQ.1) WRITE(6,999) NFIC,NLIGN
      IF(LNG.EQ.2) WRITE(6,1999) NFIC,NLIGN
999   FORMAT(1X,'UNITE LOGIQUE ',1I2,'   ERREUR LIGNE ',1I6)
1999  FORMAT(1X,'LOGICAL UNIT ',1I2,'   ERROR LINE ',1I6)
900   CONTINUE
      RETOUR = .TRUE.
!
1000  CONTINUE
!
!-----------------------------------------------------------------------
!
      RETURN
      END