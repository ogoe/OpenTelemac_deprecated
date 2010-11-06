C                       *****************
                        SUBROUTINE CG1112
C                       *****************
C
     *(X,DIM1,DIM2,IKLE,NELEM,NELMAX)
C
C***********************************************************************
C BIEF VERSION 5.1           09/12/94    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C FONCTION : CHANGEMENT DE DISCRETISATION POUR UN VECTEUR
C            ICI PASSAGE DE 11 A 12
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      X         |<-- | VECTEUR A MODIFIER.
C |      IKLE      | -->| TABLE DE CONNECTIVITE.
C |      NELEM     | -->| NOMBRE D'ELEMENTS.
C |      NELMAX    | -->| NOMBRE MAXIMUM D'ELEMENTS.
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR :
C
C  SOUS-PROGRAMME APPELE :
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER         , INTENT(IN)    :: NELEM,NELMAX,DIM1,DIM2
      DOUBLE PRECISION, INTENT(INOUT) :: X(DIM1,DIM2)
      INTEGER         , INTENT(IN)    :: IKLE(NELMAX,4)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IDIM
C
      DOUBLE PRECISION TIERS
C
C-----------------------------------------------------------------------
C
      TIERS = 1.D0/3.D0
C
C-----------------------------------------------------------------------
C
*VOCL LOOP,NOVREC
CDIR$ IVDEP
      DO 20 IDIM  = 1 , DIM2
      DO 10 IELEM = 1 , NELEM
C
        X(IKLE(IELEM,4),IDIM) = TIERS * ( X(IKLE(IELEM,1),IDIM)
     *                                  + X(IKLE(IELEM,2),IDIM)
     *                                  + X(IKLE(IELEM,3),IDIM) )
C
10    CONTINUE
20    CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
