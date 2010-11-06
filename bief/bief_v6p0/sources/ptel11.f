C                       *****************
                        SUBROUTINE PTEL11
C                       *****************
C
     *(XEL,X,IKLE,NELMAX,NELEM)
C
C***********************************************************************
C BIEF VERSION 5.1           10/01/95    J-M HERVOUET (LNH) 30 87 80 18
C
C***********************************************************************
C
C FONCTION : PASSAGE D'UN VECTEUR PAR POINTS A UN VECTEUR PAR
C            ELEMENTS.
C
C            CAS D'UN TRIANGLE P1
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      XEL       |<-- |  VECTEUR SUR LES ELEMENTS
C |      X         | -->|  VECTEUR PAR POINTS
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : NEANT
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX,IKLE(NELMAX,3)
C
C-----------------------------------------------------------------------
C
C     STRUCTURES DE VECTEURS
C
      DOUBLE PRECISION, INTENT(IN)    :: X(*)
      DOUBLE PRECISION, INTENT(INOUT) :: XEL(NELMAX,3)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM
C
C-----------------------------------------------------------------------
C
      DO 10 IELEM = 1,NELEM
C
        XEL(IELEM,1)=X(IKLE(IELEM,1))
        XEL(IELEM,2)=X(IKLE(IELEM,2))
        XEL(IELEM,3)=X(IKLE(IELEM,3))
C
10    CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
