C                       *************************
                        INTEGER FUNCTION DIM1_EXT
C                       *************************
C
     *(IELM1,IELM2,STO,TYPEXT)
C
C***********************************************************************
C BIEF VERSION 6.0       05/02/2010    J-M HERVOUET (LNH) 01 30 87 80 18
C***********************************************************************
C
C FONCTION : DONNE LA PREMIERE DIMENSION DES TERMES EXTRA-DIAGONAUX
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
      USE BIEF, EX_DIM1_EXT => DIM1_EXT
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
      INTEGER IELMX,N
C
C-----------------------------------------------------------------------
C
      IELMX = 10*(IELM1/10)
C
      IF(TYPEXT.EQ.'0') THEN
C
C        NOT 0 TO ENABLE BOUND CHECKING
         DIM1_EXT = 1
C
      ELSEIF(STO.EQ.1) THEN
C
C        CLASSICAL EBE STORAGE
C
         DIM1_EXT =NBMPTS(IELMX)
C
      ELSEIF(STO.EQ.3) THEN
C
C        EDGE-BASED STORAGE
C
         IF(TYPEXT.EQ.'S') THEN
           DIM1_EXT = NBSEG(IELM1)
         ELSE
           DIM1_EXT = NBSEG(IELM1) + NBSEG(IELM2)
           N=MAX(NBPEL(IELM1),NBPEL(IELM2))
     *      -MIN(NBPEL(IELM1),NBPEL(IELM2))
           IF(N.GE.2) THEN
C            SOME SEGMENTS LINK ONLY E.G. QUADRATIC POINTS AND
C            WILL NOT BE CONSIDERED IN A RECTANGULAR MATRIX
C            THIS IS THE CASE WITH 3 SEGMENTS IN QUADRATIC TRIANGLE
             DIM1_EXT=DIM1_EXT-N*(N-1)*NBMPTS(IELMX)/2
           ENDIF
         ENDIF
C
      ELSE
C
        IF(LNG.EQ.1) WRITE(LU,100) STO
        IF(LNG.EQ.2) WRITE(LU,101) STO
100     FORMAT(1X,'DIM1_EXT : STOCKAGE NON PREVU : ',1I6)
101     FORMAT(1X,'DIM1_EXT : UNKNOWN TYPE OF STORAGE: ',1I6)
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
