C                       *****************
                        SUBROUTINE MASKTF
C                       *****************
C
     *(MASKEL,HN,HMIN,IKLE,NELEM,NPOIN)
C
C***********************************************************************
C  BIEF VERSION 5.9     20/10/08      J-M HERVOUET (LNHE) 01 30 87 80 18
C
C***********************************************************************
C
C     FUNCTION :
C
C     MASKING DRY ELEMENTS (MASKing Tidal Flats)
C
C     ALGORITHM :
C
C     HERE SIMPLE ALGORITHM: AN ELEMENT IS DRY IF ONE OF ITS HEIGHTS IS
C                            LESS THAN HMIN
C
C     USED BY SISYPHE (IN THIS CASE THE DEPTH IS GIVEN BY TELEMAC, SO
C     A SOPHISTICATED ALGORITHM TO PRESERVE THE DYNAMICS OF FLOODING
C     IS NOT NEEDED.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   MASKEL       |<-- | TABLEAU DE MASQUAGE DES ELEMENTS             |
C |                |    |  =1. : NORMAL   =0. : ELEMENT MASQUE.        |
C |   IKLE         | -->| TABLE DE CONNECTIVITE.                       |
C |   NELEM        | -->| NOMBRE D'ELEMENTS.                           |
C |   NPOIN        | -->| NOMBRE DE NOEUDS.                            |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C***********************************************************************
C                                                                      *
C APPELE PAR: SISYPHE
C                                                                      *
C***********************************************************************
C
      IMPLICIT NONE
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NELEM,NPOIN
      INTEGER, INTENT(IN)             :: IKLE(NELEM,3)
      DOUBLE PRECISION, INTENT(IN)    :: HN(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: HMIN
      DOUBLE PRECISION, INTENT(INOUT) :: MASKEL(NELEM)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,I1,I2,I3      
C
C-----------------------------------------------------------------------
C
      DO IELEM = 1,NELEM
        I1 = IKLE(IELEM,1)
        I2 = IKLE(IELEM,2)
        I3 = IKLE(IELEM,3)
        IF(HN(I1).LE.HMIN.OR.HN(I2).LE.HMIN.OR.HN(I3).LE.HMIN) THEN
          MASKEL(IELEM) = 0.D0
        ENDIF
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
