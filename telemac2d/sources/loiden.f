C                       *****************
                        SUBROUTINE LOIDEN
C                       *****************
C
     * (YAM,YS,PHI,DEB,G)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.2     19/04/96     V. GUINOT   (LHF)
C              MODIFIE LE     03/10/96  J.-M. HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C      FONCTION:   LOI DE DEBIT D'UN SEUIL DENOYE.
C      =========
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   YAM          | -->| COTE AMONT.
C |   YS           | -->| COTE DU SEUIL.
C |   PHI          | -->| COEFFICIENT DE DEBIT DU SEUIL.
C |   DEB          |<-- | DEBIT DU SEUIL.
C |   G            | -->| GRAVITE.
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : CLSING
C
C SOUS-PROGRAMMES APPELES : NEANT
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C       
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION, INTENT(INOUT) :: DEB
      DOUBLE PRECISION, INTENT(IN)    :: G,YAM,PHI,YS
C       
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTRINSIC SQRT
C
C-----------------------------------------------------------------------
C
      IF(YAM.GT.YS) THEN
        DEB=PHI*SQRT(2.D0*G)*(YAM-YS)**1.5D0
      ELSE
        DEB=0.D0
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
