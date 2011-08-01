C                      **************************
                       SUBROUTINE MEAN_GRAIN_SIZE
C                      **************************
C
C***********************************************************************
C SISYPHE VERSION 6.0
C                                                 BUI MINH DUC NOV. 2002
C
C                                                
C COPYRIGHT EDF-BAW-IFH   
C***********************************************************************
C
C  FONCTION  : Geometric mean grain sizes of active-layer and under-layer 
C
C     SUBROUTINE A REMPLIR PAR l'UTILISATEUR
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |                |    |
C |    AVAIL       |--> | SEDIMENT FRACTION FOR EACH LAYER, CLASS, POINT
CCM-v5p3
C |    AVAIL(NPOIN,NLAYER,NSICLA)
C |    NLAYER      |--> |  NUMBER OF LAYER FOR EACH POINT
C |    NSICLA      |--> |  NUMBER OF SIZE-CLASSES OF BED MATERIAL
C |    NPOIN       |--> |  NUMBER OF NODES
CCM-v5p3
C |    ACLADM      |<-- | MEAN DIAMETER OF THE ACTIVE LAYER
C |    UNLADM      |<-- | MEAN DIAMETER OF THE ACTIVE STRATUM 
C |                       = MEAN OF ALL DIFFERENT BED MATERIAL SIZES  
C |________________|____|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C PROGRAMME APPELANT : SISYPHE (CALLED AT EACH TIME STEP)
C PROGRAMMES APPELES : 
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_SISYPHE
C
      IMPLICIT NONE
      INTEGER I , J
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C-----------------------------------------------------------------------
C 
C  UNLADM IS NEEDED FOR HUNZIKER 
C 
      DO J=1,NPOIN   
        ACLADM%R(J) = 0.D0
        UNLADM%R(J) = 0.D0
        IF(NSICLA.GT.1) THEN
         DO I=1,NSICLA
          IF(AVAIL(J,1,I).GT.0.D0) THEN
            ACLADM%R(J) = ACLADM%R(J) + FDM(I)*AVAIL(J,1,I)
            UNLADM%R(J) = UNLADM%R(J) + FDM(I)*AVAIL(J,2,I)
          ENDIF
         ENDDO
        ENDIF
        IF(ACLADM%R(J).LE.0.D0) ACLADM%R(J) = FDM(1)
        IF(UNLADM%R(J).LE.0.D0) UNLADM%R(J) = ACLADM%R(J)
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
