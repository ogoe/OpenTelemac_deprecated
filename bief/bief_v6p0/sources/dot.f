C                       *****************************
                        DOUBLE PRECISION FUNCTION DOT
C                       *****************************
C
     *(NPOIN,X,Y)
C
C***********************************************************************
C BIEF VERSION 5.1           18/08/94    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : PRODUIT SCALAIRE DES VECTEURS X ET Y DE TAILLE NPOIN
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |    NPOIN       | -->| TAILLE DE X ET Y
C |    X , Y       | -->| TABLEAUX CONTENANT LES VECTEURS
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NPOIN
      DOUBLE PRECISION, INTENT(IN) :: X(NPOIN),Y(NPOIN)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I
C
C-----------------------------------------------------------------------
C
      DOT = 0.D0
C
      DO I = 1 , NPOIN
       DOT = DOT + X(I) * Y(I)
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
