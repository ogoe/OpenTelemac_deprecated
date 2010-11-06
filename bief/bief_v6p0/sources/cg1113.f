C                       *****************
                        SUBROUTINE CG1113
C                       *****************
C
     *(X,DIM1,DIM2,IKLE,NELEM,NELMAX)
C
C***********************************************************************
C BIEF VERSION 5.9         06/02/08    J-M HERVOUET (LNH) 01 30 87 80 18
C***********************************************************************
C
C FONCTION : CHANGEMENT DE DISCRETISATION POUR UN VECTEUR
C            ICI PASSAGE DE 11 A 13 (LINEAIRE A QUADRATIQUE)
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
      INTEGER         , INTENT(IN)    :: IKLE(NELMAX,6)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IDIM
C
C-----------------------------------------------------------------------
C
*VOCL LOOP,NOVREC
CDIR$ IVDEP
      DO IDIM  = 1 , DIM2
      DO IELEM = 1 , NELEM
C
        X(IKLE(IELEM,4),IDIM) = 0.5D0 * ( X(IKLE(IELEM,1),IDIM)
     *                                  + X(IKLE(IELEM,2),IDIM) )
        X(IKLE(IELEM,5),IDIM) = 0.5D0 * ( X(IKLE(IELEM,2),IDIM)
     *                                  + X(IKLE(IELEM,3),IDIM) )
        X(IKLE(IELEM,6),IDIM) = 0.5D0 * ( X(IKLE(IELEM,3),IDIM)
     *                                  + X(IKLE(IELEM,1),IDIM) )
C
      ENDDO
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
