C                       ***********************
                        LOGICAL FUNCTION INCLUS
C                       ***********************
C
     *( C1 , C2 )
C
C***********************************************************************
C BIEF VERSION 5.1        17/08/94   J.M. HERVOUET (LNH)   30 87 80 18
C
C***********************************************************************
C
C FONCTION  : CONTROLE QU'UNE CHAINE DE CARACTERES EST INCLUSE DANS
C             UNE AUTRE :
C
C             INCLUS=.TRUE. EST EQUIVALENT A 'C2 INCLUSE DANS C1'
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   C2           | -->| CHAINE CENSEE CONTENIR C1
C |   C1           | -->| CHAINE RECHERCHEE DANS C2
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
      INTEGER I,LC1,LC2
C
      INTRINSIC LEN
C
C-----------------------------------------------------------------------
C
      INCLUS = .FALSE.
C
      LC1 = LEN(C1)
      LC2 = LEN(C2)
      IF(LC2.GT.LC1) GO TO 1000
C
      I = 0
10    I = I + 1
      IF(I.GT.LC1-LC2+1) GO TO 1000
C
      IF(C1(I:I+LC2-1).NE.C2(1:LC2)) GO TO 10
C
      INCLUS = .TRUE.
C
1000  CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
