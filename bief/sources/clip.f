C                       ***************
                        SUBROUTINE CLIP
C                       ***************
C
     *(F,XMIN,CLPMIN,XMAX,CLPMAX,NPOIN)
C
C***********************************************************************
C BIEF VERSION 5.1           12/01/95    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : LIMITATION DES VALEURS CONTENUES DANS LE TABLEAU F
C
C            MINIMUM DE F : XMIN (SI CLPMIN=.TRUE.)
C            MAXIMUM DE F : XMAX (SI CLPMAX=.TRUE.)
C
C
C IMPORTANTE NOTE : SI NPOIN EST NEGATIF, ON TRAITERA - NPOIN VALEURS
C                   SI NPOIN EST POSITIF ON PRENDRA LA TAILLE DE F.
C
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |    F           |<-->| TABLEAU DES VALEURS
C |    XMIN        | -->| VALEUR MIN
C |    CLPMIN      | -->| LOGIQUE QUI DECIDE DU CLIPPING
C |    XMAX        | -->| VALEUR MAX
C |    CLPMAX      | -->| LOGIQUE QUI DECIDE DU CLIPPING
C |    NPOIN       | -->| NOMBRE TOTAL DE POINTS DU DOMAINE
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C***********************************************************************
C
      USE BIEF, EX_CLIP => CLIP
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ)  , INTENT(INOUT) :: F
      DOUBLE PRECISION, INTENT(IN)    :: XMIN,XMAX
      LOGICAL         , INTENT(IN)    :: CLPMIN,CLPMAX
      INTEGER         , INTENT(IN)    :: NPOIN
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER NP
C
C-----------------------------------------------------------------------
C
      IF(F%TYPE.EQ.2) THEN
C       F EST UNE STRUCTURE DE VECTEUR
        IF(NPOIN.LT.0) THEN
          NP = - NPOIN
        ELSE
          NP = F%DIM1
        ENDIF
        IF(CLPMIN) CALL OV('X=+(Y,C)',F%R,F%R, F%R , XMIN , NP )
        IF(CLPMAX) CALL OV('X=-(Y,C)',F%R,F%R, F%R , XMAX , NP )
      ELSE
        IF(LNG.EQ.1) WRITE(LU,*) F%NAME,' N''EST PAS UN VECTEUR'
        IF(LNG.EQ.2) WRITE(LU,*) F%NAME,' IS NOT A VECTOR'
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
