CV 04/05/2010  modif pour Fredsoe: les conc d'equilibre doivent-être multipliees par
C avai: sinon le calcul des concentrations en entree du domaine est faux (ie ne tient 
C pas compte des avai)
C
      ! ****************************** !
        SUBROUTINE SUSPENSION_EROSION  !
      ! ****************************** !
C   
     *(TAUP,HN,ACLADM,AVA,NPOIN,CHARR,XMVE,XMVS,GRAV,HMIN,XWC,
     * ZERO,ZREF,AC,FLUER,CSTAEQ,QSC,ICQ,DEBUG)
C
C**********************************************************************
C SISYPHE VERSION 6.0  17/09/09               J-M HERVOUET + C VILLARET            
C
C**********************************************************************
C
C          ! ================================================= !
C          ! Computation of the flux of deposition and erosion !
C          ! ================================================= !
C
C
C COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT
C**********************************************************************C
C                                                                      C
C                 SSSS I   SSSS Y   Y PPPP  H   H EEEEE                C
C                S     I  S      Y Y  P   P H   H E                    C
C                 SSS  I   SSS    Y   PPPP  HHHHH EEEE                 C
C                    S I      S   Y   P     H   H E                    C
C                SSSS  I  SSSS    Y   P     H   H EEEEE                C
C                                                                      C
C----------------------------------------------------------------------C
C                             ARGUMENTS                                C
C .________________.____.______________________________________________C
C |      NOM       |MODE|                   ROLE                       C
C |________________|____|______________________________________________C
C |    TOB         | => |  FLUX DE DEPOT                               C
C |    CF          | => |  
C |    HN          | => |  HAUTEUR D'EAU
C |________________|____|______________________________________________C
C                    <=  Can't be change by the user                   C
C                    =>  Can be changed by the user                    C 
C ---------------------------------------------------------------------C
!                                                                      !
! CALLED BY SUSPENSION_COMPUTATION                                     !
!                                                                      !
! CALL     SUSPENSION_FREDSOE                                          !
!          SUSPENSION_BIJKER                                           !
!                                                                      !
!======================================================================!
!======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS               !
!======================================================================!
!======================================================================!

      ! 1/ MODULES
      ! ----------
      USE INTERFACE_SISYPHE, EX_SUSPENSION_EROSION=>SUSPENSION_EROSION
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU         

      ! 2/ GLOBAL VARIABLES
      ! -------------------
      TYPE (BIEF_OBJ),  INTENT(IN)    :: TAUP,HN,ACLADM,ZREF,QSC
      INTEGER,          INTENT(IN)    :: NPOIN,DEBUG,ICQ
      LOGICAL,          INTENT(IN)    :: CHARR
      DOUBLE PRECISION, INTENT(IN)    :: XMVE,XMVS,GRAV,HMIN,XWC,ZERO
      DOUBLE PRECISION, INTENT(IN)    :: AVA(NPOIN)
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: FLUER,CSTAEQ
      DOUBLE PRECISION, INTENT(INOUT) :: AC
      
      ! 3/ LOCAL VARIABLES
      ! -------------------

      INTEGER I

!
!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!
! 
!
!  COMPUTATION OF THE NEAR BED EQUILIBRIUM CONC --> CSTAEQ (mean diameter)
!           
      IF(ICQ.EQ.1) THEN 
C
        IF(DEBUG > 0) WRITE(LU,*) 'SUSPENSION_FREDSOE'
        CALL SUSPENSION_FREDSOE(ACLADM,TAUP,NPOIN,
     &                          GRAV,XMVE,XMVS,ZERO,AC,CSTAEQ)
        IF(DEBUG > 0) WRITE(LU,*) 'END SUSPENSION_FREDSOE'
C
C       CALL OS('X=CYZ   ', X=FLUER, Y=CSTAEQ, Z=AVA, C=XWC)
C 04/05/2010
C debut modifs CV ...
CV        DO I=1,NPOIN
CV          FLUER%R(I)=XWC*CSTAEQ%R(I)*AVA(I)
CV        ENDDO
          DO I=1,NPOIN
            CSTAEQ%R(I)=CSTAEQ%R(I)*AVA(I)
          ENDDO
          CALL OS('X=CY    ', X=FLUER, Y=CSTAEQ, C=XWC)            
C ... fin modifs CV
C  
      ELSEIF(ICQ.EQ.2) THEN
C
        IF(DEBUG > 0) WRITE(LU,*) 'SUSPENSION_BIJKER'
        CALL SUSPENSION_BIJKER(TAUP,HN,NPOIN,CHARR,QSC,ZREF,
     &                         ZERO,HMIN,CSTAEQ,XMVE)
        IF(DEBUG > 0) WRITE(LU,*) 'END SUSPENSION_BIJKER'
C       pas de  multiplication par ava car ava est déjà pris en compte 
C       dans le taux de transport par charriage
        CALL OS('X=CY    ', X=FLUER, Y=CSTAEQ, C=XWC)            
C
      ENDIF    
!
!======================================================================!
!======================================================================!
!
      RETURN      
      END
