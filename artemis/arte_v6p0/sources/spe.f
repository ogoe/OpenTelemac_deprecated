C                       *****************************
                        DOUBLE PRECISION FUNCTION SPE
C                       *****************************
C
     *(F)
C
C***********************************************************************
C
C  ARTEMIS    VERSION 2.0  01/06/93   F. LEPEINTRE (LNH) 01 30 87 78 54
C             VERSION 5.1  04/06/99   D. AELBRECHT (LNH) 01 30 87 74 12
C
C***********************************************************************
C
C      FONCTION:    CALCULE LA DENSITE D'ENERGIE SUIVANT LA FORMULE
C                   DE GODA:RANDOM SEA AND DESIGN OF MARITIME STRUCTURES
C                           UNIVERSITY OF TOKYO PRESS - 1985
C
C
C SPE(F) = BETA * HS**2 * FP**4/F**5 * EXP(-1.25*(FP/F)**4) *
C          GAM**(EXP(-0.5*( (F-FP)/(SIGMA*FP) )**2)
C
C
C AVEC BETA = 0.0624 / ( 0.230+O.0336*GAM-0.185 / (1.9+GAM) )
C      HS   = 1 (HAUTEUR SIGNIFICATIVE)
C      FP : FREQUENCE DE PIC
C      GAM : VALEUR CARACTERISTIQUE DU SPECTRE
C              (GAM=1 PIERSON-MOSKOWITZ, GAM=3.3 JONSWAP)
C      F : FREQUENCE OU ON CALCULE LA DENSITE
C      SIGMA = 0.07 SI F < FP  OU F = FP
C      SIGMA = 0.09 SI F > FP
C
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   F            | -->|  FREQUENCE OU L'ON CALCULE LA DENSITE        |
C |                |    |  D'ENERGIE                                   |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : PERALE
C
C***********************************************************************
C
C       USE INTERFACE_ARTEMIS 
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      DOUBLE PRECISION F,SIGMA
C
      DOUBLE PRECISION FP,GAM,DELTA
      COMMON /COEFHE/ FP,GAM,DELTA
C
      INTRINSIC EXP
C
C-----------------------------------------------------------------------
C
      IF (F.LE.FP) THEN
        SIGMA = 0.07D0
      ELSE
        SIGMA = 0.09D0
      ENDIF
C
CC    DELTA = 0.0624D0 * FP**4 /
CC   *       ( 0.230D0 + 0.0336D0*GAM - 0.185D0 / (1.9D0+GAM) )
C    
C     DELTA EST CALCULEE DANS PERALE
C
      IF ( F.GE.1.D-4*FP) THEN
         SPE = DELTA/F**5 * EXP(-1.25D0*(FP/F)**4) *
     *         GAM** ( EXP( -0.5D0*( ( (F-FP)/(SIGMA*FP) ) **2 ) ) )
      ELSE
         SPE = 0.D0
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
