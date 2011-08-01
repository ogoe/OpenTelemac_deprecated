C                       *****************
                        SUBROUTINE CNTPRE
C                       *****************
C
     *(DAM,NPOIN,IPRECO,IPREC2)
C
C***********************************************************************
C
C  ARTEMIS VERSION 5.1  02/06/99   D. AELBRECHT (LNH) 01 30 87 74 12 
C
C  LINKED TO BIEF VERS. 5.0          J-M HERVOUET (LNH) 01 30 87 80 18
C
C***********************************************************************
C
C      FONCTION:    INHIBE LE PRECONDITIONNEMENT DIAGONAL SI UN DES
C                   ELEMENTS DE DAM EST NEGATIF OU NUL.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   DAM          | -->|  DIAGONALE DE LA MATRICE                     |
C |   NPOIN        | -->|  NOMBRE DE POINTS DU MAILLAGE                |
C |   IPRECO       | -->|  PRECONDITIONNEMENT DEMANDE PAR L'UTILISATEUR|
C |   IPREC2       |<-- |  PRECONDITIONNEMENT UTILISE                  |
C |________________|____|______________________________________________|
C
C-----------------------------------------------------------------------
C  VALEUR DE IPREC2   I                  SIGNIFICATION
C  OU DE IPRECO       I
C-----------------------------------------------------------------------
C        1            I  RIEN
C        2            I  PRECONDITIONNEMENT DIAGONAL AVEC LA DIAGONALE
C                     I  DE LA MATRICE.
C        3            I  PRECONDITIONNEMENT DIAGONAL AVEC LA MATRICE
C                     I  CONDENSEE.
C        5            I  AUTRE PRECONDITIONNEMENT (NON ATTRIBUE)
C                     I
C        7            I  PRECONDITIONNEMENT DE CROUT PAR ELEMENT
C                     I  (NON PROGRAMME)
C
C-----------------------------------------------------------------------
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : BERKHO
C
C***********************************************************************
C              
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INTEGER NPOIN,IPRECO,IPREC2,I
C
      DOUBLE PRECISION DAM(NPOIN)
C
C-----------------------------------------------------------------------
C
      IF (IPRECO.EQ.0) IPRECO = 1
      IPREC2 = IPRECO
C
      IF (MOD(IPRECO,2).EQ.0.OR.MOD(IPRECO,3).EQ.0) THEN
C
         DO 10 I=1,NPOIN
            IF (DAM(I).LE.0.D0) THEN
20             CONTINUE
               IF(MOD(IPREC2,2).EQ.0) THEN
                 IPREC2 = IPREC2/2
                 GO TO 20
               ENDIF
21             CONTINUE
               IF(MOD(IPREC2,3).EQ.0) THEN
                 IPREC2 = IPREC2/3
                 GO TO 21
               ENDIF
               IF (LNG.EQ.1) WRITE(LU,100)
               IF (LNG.EQ.2) WRITE(LU,101)
100     FORMAT(1X,'CNTPRE (ARTEMIS) : PRECONDITIONNEMENT DIAGONAL NON AP
     *PLIQUE (UN ELEMENT DIAGONAL DE LA MATRICE EST NEGATIF OU NUL)')
101     FORMAT(1X,'CNTPRE (ARTEMIS) : DIAGONAL SCALING NOT APPLIED (ONE
     *COEFFICIENT OF THE MATRIX DIAGONAL IS NEGATIVE OR ZERO)')
               GOTO 30
            ENDIF
10       CONTINUE
C
      ELSEIF (MOD(IPRECO,5).EQ.0) THEN
C
         DO 40 I=1,NPOIN
            IF (ABS(DAM(I)).LE.1.D-6) THEN
50             CONTINUE
               IF(MOD(IPREC2,5).EQ.0) THEN
                 IPREC2 = IPREC2/5
                 GO TO 50
               ENDIF
               IF (LNG.EQ.1) WRITE(LU,200)
               IF (LNG.EQ.2) WRITE(LU,201)
200     FORMAT(1X,'CNTPRE (ARTEMIS) : PRECONDITIONNEMENT DIAGONAL NON AP
     *PLIQUE (UN ELEMENT DIAGONAL DE LA MATRICE EST NUL)')
201     FORMAT(1X,'CNTPRE (ARTEMIS) : DIAGONAL SCALING NOT APPLIED (ONE
     *COEFFICIENT OF THE MATRIX DIAGONAL IS ZERO)')
               GOTO 30
            ENDIF
40       CONTINUE
C
      ENDIF
C
30    CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
