C                       *******************************
                        DOUBLE PRECISION FUNCTION P_DOT
C                       *******************************
C
     *(NPOIN,X,Y,FAC)
C
C***********************************************************************
C BIEF VERSION 5.5           24/04/97    J-M HERVOUET (LNH) 30 87 80 18
C COPYRIGHT 1997              AFTER REINHARD HINKELMANN (HANNOVER UNI.)
C***********************************************************************
C
C FONCTION : PRODUIT SCALAIRE DES VECTEURS X ET Y DE TAILLE NPOIN
C            AVEC PRISE EN COMPTE DU PARALLELISME.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |    NPOIN       | -->| TAILLE DE X ET Y
C |    X , Y       | -->| TABLEAUX CONTENANT LES VECTEURS
C |    FAC         | -->| FAC=1/(nombre de domaines voisins du point)
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C***********************************************************************
C
      USE BIEF, EX_P_DOT => P_DOT
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      INTEGER, INTENT(IN) :: NPOIN
C
      DOUBLE PRECISION, INTENT(IN) :: X(NPOIN),Y(NPOIN),FAC(NPOIN)
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      INTEGER I
C
C-----------------------------------------------------------------------
C
      P_DOT = 0.D0
C
      DO 10 I = 1 , NPOIN
       P_DOT = P_DOT + X(I) * Y(I) * FAC(I)
10    CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
