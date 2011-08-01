C                       ***********************
                        LOGICAL FUNCTION INCLU2
C                       ***********************
C
     *( C1 , C2 )
C
C***********************************************************************
C BIEF VERSION 5.1        17/08/94   J.M. HERVOUET (LNH)   30 87 80 18
C
C***********************************************************************
C
C FONCTION  : CONTROLE QU'UN MOT EST INCLUS DANS UNE LISTE DE MOTS :
C
C          INCLU2=.TRUE. EST EQUIVALENT A 'MOT C2 INCLUS DANS LISTE C1'
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   C1           | -->| LISTE DE MOTS SEPARES PAR AUTRE CHOSE QUE    |
C |                |    |   A-Z ET 0-9
C |   C2           | -->| MOT RECHERCHE DANS C1                        |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      CHARACTER*(*) C1 , C2
C
      INTEGER I,IC1,LC1,LC2,IMAX
C
      LOGICAL FLAG
C
      INTRINSIC LEN
C
C-----------------------------------------------------------------------
C
      INCLU2 = .FALSE.
C
      LC1 = LEN(C1)
      LC2 = LEN(C2)
      IMAX = LC1-LC2
C
      IF(IMAX.GE.0) THEN
C
         DO 10 I = 0,IMAX
            IF(C1(I+1:I+LC2).EQ.C2(1:LC2)) THEN
               FLAG = .TRUE.
               IF (I.NE.0) THEN
                  IC1 = ICHAR(C1(I:I))
                  IF ((IC1.GE.48.AND.IC1.LE.57).OR.
     *                (IC1.GE.65.AND.IC1.LE.90)) FLAG = .FALSE.
               ENDIF
               IF (I.NE.IMAX) THEN
                  IC1 = ICHAR(C1(I+LC2+1:I+LC2+1))
                  IF ((IC1.GE.48.AND.IC1.LE.57).OR.
     *                (IC1.GE.65.AND.IC1.LE.90)) FLAG = .FALSE.
               ENDIF
               INCLU2 = INCLU2.OR.FLAG
            ENDIF
10       CONTINUE
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
