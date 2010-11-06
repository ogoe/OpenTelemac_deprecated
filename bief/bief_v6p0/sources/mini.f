C                       ***************
                        SUBROUTINE MINI
C                       ***************
C
     *( XMIN , IMIN , X , NPOIN )
C
C***********************************************************************
C BIEF VERSION 5.1       17/08/94    E. PELTIER   (LNH)
C***********************************************************************
C
C  FONCTION : RECHERCHE DE LA PLUS PETITE VALEUR DANS UN TABLEAU X
C             DE DIMENSION NPOIN
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |      XMIN      |<-- | MINIMUM TROUVE
C |      IMIN      |<-- | INDICE DU MINIMUM
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
      INTEGER, INTENT(INOUT)          :: IMIN
      DOUBLE PRECISION, INTENT(INOUT) :: XMIN
      DOUBLE PRECISION, INTENT(IN)    :: X(NPOIN)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I
C
C-----------------------------------------------------------------------
C
      XMIN = X(1)
      IMIN = 1
C
      DO 10 I = 2, NPOIN
C
        IF(X(I).LT.XMIN) THEN
          IMIN = I
          XMIN = X(I)
        ENDIF
C
10    CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
