C                       *********************************
                        SUBROUTINE CALC_CRIT_BOTTOM_ANGLE
C                       *********************************
     *( BETA_CRIT )
C
C
C***********************************************************************
C SISYPHE VERSION 6.0                   16/02/11       O. GOETHEL
C
C***********************************************************************
C
C     FUNCTION  : Calculate the critical bottom angle
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   BETA_CRIT    |<-- | Critical Bottom Angle
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C***********************************************************************
C
C
      USE DECLARATIONS_SISYPHE, ONLY  : NPOIN,ACLADM
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C
C
      INTEGER        :: IPOIN
      DOUBLE PRECISION, POINTER :: DM
      DOUBLE PRECISION :: BETA_CRIT(NPOIN)
      DOUBLE PRECISION :: PI,FA
C
C
      PI = ACOS(-1.D0)
      BETA_CRIT=0.D0
C
      CALL MEAN_GRAIN_SIZE()
C
      DO IPOIN=1,NPOIN
C
       DM => ACLADM%R(IPOIN)      
C
       IF(DM.LT.0.00006 .AND. DM.GE.0.000002)  THEN
       
        FA = (DM-0.000002D0)/0.000058D0 * 5.D0 + 25.D0

       ELSEIF(DM.LT.0.005 .AND. DM.GE.0.00006)THEN
       
        FA = (DM-0.00006D0)/0.00494D0 * 2.D0 + 30.D0
       
       ELSEIF(DM.LT.0.01 .AND. DM.GE.0.005) THEN
       
        FA = (DM-0.005D0)/0.005D0 * 3.D0 + 32.D0
       
       ELSEIF(DM.LT.0.05 .AND. DM.GE.0.01)  THEN
       
        FA = (DM-0.01D0)/0.04D0 * 2.D0 + 35.D0
       
       ELSEIF(DM.LE.0.1 .AND. DM.GE.0.05)   THEN
       
        FA = (DM-0.05D0)/0.05D0 * 3.D0 + 37.D0
       
       ELSE
       
       write(LU,*) 'CRITICAL BOTTOM ANGLE CAN NOT BE DETERMINED'
       STOP
       
       ENDIF
C
       FA = FA * PI/180.D0
C
       BETA_CRIT(IPOIN) = FA
C
      ENDDO
C
      END SUBROUTINE CALC_CRIT_BOTTOM_ANGLE