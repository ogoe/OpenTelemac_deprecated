C                       ***************
                        SUBROUTINE MAXI
C                       ***************
C
     *( XMAX , IMAX , X , NPOIN )
C
C***********************************************************************
C BIEF VERSION 5.1       17/08/94    E. PELTIER   (LNH)
C***********************************************************************
C
C  FONCTION : RECHERCHE DE LA PLUS GRANDE VALEUR DANS UN TABLEAU X
C             DE DIMENSION NPOIN
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |      XMAX      |<-- | MAXIMUM TROUVE
C |      IMAX      |<-- | INDICE DU MAXIMUM
C |      X         | -->| TABLEAU DES VALEURS
C |      NPOIN     | -->| DIMENSION DU TABLEAU
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C PROGRAMME APPELANT :
C PROGRAMMES APPELES : AUCUN
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NPOIN 
      INTEGER, INTENT(INOUT)          :: IMAX 
      DOUBLE PRECISION, INTENT(INOUT) :: XMAX
      DOUBLE PRECISION, INTENT(IN)    :: X(NPOIN)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I     
C
C-----------------------------------------------------------------------
C
      XMAX = X(1)
      IMAX = 1
C
      DO 10 I = 2 , NPOIN
C
        IF(X(I).GT.XMAX) THEN
          IMAX = I
          XMAX = X(I)
        ENDIF
C
10    CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
