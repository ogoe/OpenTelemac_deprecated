C                       *******************************
                        DOUBLE PRECISION FUNCTION SOMME
C                       *******************************
C
     *( X , NPX )
C
C***********************************************************************
C BIEF VERSION 5.1           08/12/94    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C   FONCTION : SOMME DES COMPOSANTES D'UN VECTEUR
C              (VOIR AUSSI SOMME2 ET SUM)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|_______________________________________________
C |      X         | -->| TABLEAU FORTRAN
C |      NPX       | -->| NOMBRE DE VALEURS A ADDITIONNER
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES   :
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INTEGER NPX,I
C
      DOUBLE PRECISION X(*)
C
C-----------------------------------------------------------------------
C
        SOMME = 0.D0
        DO 10 I = 1 , NPX
          SOMME = SOMME + X(I)
10      CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
