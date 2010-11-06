C                       *************************
                        INTEGER FUNCTION DIM2_EXT
C                       *************************
C
     *(IELM1,IELM2,STO,TYPEXT)
C
C***********************************************************************
C BIEF VERSION 5.1           12/12/94    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C FONCTION : DONNE LA DEUXIEME DIMENSION DES TERMES EXTRA-DIAGONAUX
C            D'UNE MATRICE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |  IELM1         | -->| TYPE DE L'ELEMENT DE LIGNE
C |  IELM2         | -->| TYPE DE L'ELEMENT DE COLONNE
C |  STO           | -->| TYPE DE STOCKAGE
C |  TYPEXT        | -->| TYPE DES TERMES EXTRA-DIAGONAUX
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C SOUS-PROGRAMME APPELE : NEANT
C
C***********************************************************************
C
      USE BIEF, EX_DIM2_EXT => DIM2_EXT
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER         , INTENT(IN) :: IELM1,IELM2,STO
      CHARACTER(LEN=1), INTENT(IN) :: TYPEXT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER NDIM
C
C-----------------------------------------------------------------------
C
      IF(TYPEXT.EQ.'0') THEN
C
         DIM2_EXT = 1
C
      ELSEIF(STO.EQ.1) THEN
C
C        NDIM IS HERE THE SECOND DIMENSION OF VECTOR
         NDIM =     NBPEL(IELM1)*NBPEL(IELM2)
     *        - MIN(NBPEL(IELM1),NBPEL(IELM2))
C        SYMMETRICAL MATRIX
         IF(IELM1.EQ.IELM2.AND.TYPEXT.EQ.'S') THEN
           NDIM = NDIM / 2
         ENDIF
         DIM2_EXT = NDIM
C
      ELSEIF(STO.EQ.3) THEN
C
C        ASSEMBLED EBE STORAGE OR EDGE BASED
         DIM2_EXT = 1
C
      ELSE
C
        IF(LNG.EQ.1) WRITE(LU,100) STO
        IF(LNG.EQ.2) WRITE(LU,101) STO
100     FORMAT(1X,'DIM2_EXT : STOCKAGE NON PREVU : ',1I6)
101     FORMAT(1X,'DIM2_EXT : UNKNOWN TYPE OF STORAGE: ',1I6)
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
    
