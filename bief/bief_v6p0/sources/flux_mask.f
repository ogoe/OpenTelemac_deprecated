C                       ********************
                        SUBROUTINE FLUX_MASK
C                       ********************
C
     *(FXMAT,NSEG,GLOSEG,SIZGLO,MASKPT)
C
C***********************************************************************
C BIEF VERSION 5.9       19/06/08     J-M HERVOUET (LNHE) 01 30 87 80 18
C
C***********************************************************************
C
C  FONCTION  : MASQUAGE DES FLUX PAR SEGMENT AVEC LES MASQUES DES
C              EXTREMITES DU SEGMENT
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|_______________________________________________
C |    FXMAT       |<-- | MATRICE DE STOCKAGE DES FLUX.
C |    NSEG        | -->| NOMBRE DE SEGMENTS DANS LE MAILLAGE.
C |    NELEM       | -->| NOMBRE D'ELEMENTS DANS LE MAILLAGE.
C |    ELTSEG      | -->| SEGMENTS OF EVERY TRIANGLE.
C |    ORISEG      | -->| ORIENTATION OF SEGMENTS OF EVERY TRIANGLE.
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT : CVDFTR
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF, EX_FLUX_MASK => FLUX_MASK
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NSEG,SIZGLO
      INTEGER, INTENT(IN)             :: GLOSEG(SIZGLO,2)
      DOUBLE PRECISION, INTENT(INOUT) :: FXMAT(NSEG)
      DOUBLE PRECISION, INTENT(IN)    :: MASKPT(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I
C
C-----------------------------------------------------------------------
C
      DO I = 1,NSEG
        FXMAT(I) = FXMAT(I) * MASKPT(GLOSEG(I,1)) * MASKPT(GLOSEG(I,2))
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
