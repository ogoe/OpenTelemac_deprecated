C                       *************************
                        SUBROUTINE CONDIM_SISYPHE
C                       *************************
C
     * (U      , V       , QU    , QV   , H    , ZF , Z ,
     *  ESOMT  , THETAWR ,  Q    , HWR  , TWR  , 
     *  X      , Y       , NPOIN , AT   , PMAREE)
C
C ATTENTION MODIF CONDIM MAI 2006
C LES VARIABLES U et V,  H, doivent être  définies 
C               ZF
C LES AUTRES SONT OPTIONNELLES 
C***********************************************************************
C SISYPHE VERSION 5.3                             E. PELTIER    11/09/95
C                                                 C. LENORMANT
C                                                 J.-M. HERVOUET
C                                                
C COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT      
C***********************************************************************
C
C     FONCTION  : VALEURS IMPOSEES
C                         - DU DEBIT VECTORIEL          QU, QV
C                         - DE LA HAUTEUR D'EAU         H
C                         - DE LA COTE DU FOND          ZF
C                         - DE LA SURFACE LIBRE         Z
C                         - DE L'EVOLUTION TOTALE       ESOMT
C                         - DU DEBIT                    Q
C                         - DE LA HAUTEUR DE HOULE      HW
C                         - DE LA PERIODE DE HOULE      TW
C                         - DE L'ANGLE D'ATTAQUE DE LA HOULE (PAR RAPPORT A L'AXE OY)
C                                                       THETAW
C
C      FUNCTION :       IMPOSED VALUES OF
C                         -  DEPTH-AVERAGED FLOW RATE  (X,Y) 
C                                              QU, QV
C                         -   WATER DEPTH               H
C                         -   BOTTOM ELEVATION          ZF
C                         -   FREE SURFACE              Z
C                         -   TOTAL BED VOLUTION        ESOMT
C                         -   FLOW RATE                 Q
C                         -   WAVE HEIGHT               HWR
C                         -   WAVE PERIOD               TWR
C                         -   WAVE DIRECTION      (with axis Oy)
C                                                       THETAWR      
C 
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   U , V        |<-- | FLOW VELOCITY COORDINATES
C |   QU , QV      |<-- | FLOW RATE COORDINATES
C |   H            |<-->| WATER DEPTH
C |   ZF           |<-->| BED ELEVATION
C |   Z            |<-->| FREE SURFACE
C |   ESOMT        |<-->| TOTAL BED EVOLUTION
C |   Q            |<-->| FLOW RATE
C |   HW           |<-->| WAVE HEIGHT (M)
C |   TW           |<-->| WAVE PERIOD (S)
C |   TETHAW       |<-->| WAVE ANGLE (DEG)
C |   X,Y          | -->| COORDINATES
C |   NPOIN        | -->| NUMBER OF 2D POINTS
C |   AT           | -->| TIME
C |   PMAREE       | -->| TIDAL PERIOD
C |________________|____|______________________________________________
C MODE : -->(INPUT ), <--(RESULT), <-->(MODIFIED VARIABLE)
C-----------------------------------------------------------------------
C CALLED BY SUBROUTINE  : SISYPHE
C 
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_SISYPHE, ONLY: HW,TW,THETAW
C
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)::NPOIN 
C
      DOUBLE PRECISION, INTENT(IN):: X(NPOIN),Y(NPOIN)
      DOUBLE PRECISION, INTENT(IN):: AT , PMAREE 
C sediment 
      DOUBLE PRECISION, INTENT(INOUT) ::  ZF(NPOIN)
      DOUBLE PRECISION, INTENT (INOUT)::  ESOMT(NPOIN)
C hydrodynamics
      DOUBLE PRECISION, INTENT(INOUT):: Z(NPOIN) , H(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT):: U(NPOIN) , V(NPOIN)
      DOUBLE PRECISION, INTENT (INOUT)::QU(NPOIN), QV(NPOIN), Q(NPOIN)
C waves 
      DOUBLE PRECISION, INTENT (INOUT):: HWR(NPOIN) , TWR(NPOIN)
      DOUBLE PRECISION, INTENT (INOUT):: THETAWR(NPOIN)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C      INTEGER I
C-----------------------------------------------------------------------
C
C  INITIALIZATION OF VARIABLES NOT READ IN RESULTS FILE
C   (REPLACEMENT OF VALUES READ IN RESULTS FILE)
C
C     ------------------------
C     THE USER SHOULD BE AWARE 
C     ++++++++++++++++++++++++
C
C     SUBROUTINE CONDIM_SISYPHE IS CALLED AT EACH TIME STEP 
C     IN ORDER TO IMPOSE A VARIABLE FORCING
C     (TIDAL CURRENT, FOR EXAMPLE)
C    
C     IT IS NOT SUFFICIENT TO PRESCRIBE THE FLOW RATE
C     THE MAIN VARIABLES ARE NOW THE 2D FLOW VELOCITY FIELD
C     AND THE FLOW DEPTH
C
C-----------------------------------------------------------------------
C  
C     WAVES, EXAMPLE WITH NO WAVES:
C
C     AMPLITUDE = 0
!     CALL OS('X=0     ',X=HW)
C     PERIOD = 1 S
!     CALL OS('X=C     ',X=TW,C=1.D0)
C     ANGLE = 0
!     CALL OS('X=0     ',X=THETAW)
C
C     AFTER SETTING HWR, TWR AND THETAWR, PLEASE ADD:
C
!     HW%TYPR    ='Q'
!     TW%TYPR    ='Q'
!     THETAW%TYPR='Q'
C
C     TO ENABLE THE CONTROL OF WAVE DATA
C
C-----------------------------------------------------------------------
C                                                
      RETURN
      END SUBROUTINE CONDIM_SISYPHE
