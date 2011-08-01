      ! ****************************** !
        SUBROUTINE SUSPENSION_EROSION_COH  
      ! ****************************** !

     &(TAUP,NPOIN,XMVE,XMVS,GRAV,VITCE,
     & PARTHENIADES,ZERO,DEBUG, 
     & FLUER, ES, TOCE_VASE, NCOUCH_TASS, DT, MS_VASE,TASS)

C**********************************************************************
C SISYPHE VERSION 6.0  31/07/08                             C. VILLARET            
C
C**********************************************************************
C
C          ! ================================================= !
C          ! Computation of the flux of deposition and erosion !
C          ! ================================================= !
C
CV nouvelle subroutine pour calcul du fluerosion en tenant compte de 
CV la structure verticale
C
C 
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
!                                                                      !                                      !
!                                                                      !
!======================================================================!
!======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS               !
!======================================================================!
!======================================================================!

      ! 1/ MODULES
      ! ----------
      USE INTERFACE_SISYPHE, EX_SUSPENSION_EROSION_COH=>
     *                          SUSPENSION_EROSION_COH
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU         

      ! 2/ GLOBAL VARIABLES
      ! -------------------
      TYPE (BIEF_OBJ),  INTENT(IN)    :: TAUP
      INTEGER,          INTENT(IN)    :: NPOIN,DEBUG
      DOUBLE PRECISION, INTENT(IN)    :: XMVE,XMVS,GRAV
      DOUBLE PRECISION, INTENT(IN)    :: VITCE
      DOUBLE PRECISION, INTENT(IN)    :: ZERO,PARTHENIADES
! pour le tassement
      DOUBLE PRECISION,  INTENT(INOUT) :: MS_VASE(NPOIN,10)
      DOUBLE PRECISION, INTENT(IN)     :: TOCE_VASE(10), DT
      INTEGER,          INTENT(IN)    :: NCOUCH_TASS
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: FLUER
      DOUBLE PRECISION, INTENT(INOUT) :: ES(NPOIN,10)
!
      LOGICAL, INTENT(IN) :: TASS

      ! 3/ LOCAL VARIABLES
      ! ------------------
      INTEGER :: I, J

      DOUBLE PRECISION :: USTARP,AUX
      DOUBLE PRECISION :: FLUER_LOC(10), QER_VASE,TEMPS, QE_COUCHE
      
!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!

      ! *************************************************  !
      ! Ia - FORMULATION FOR COHESIVE SEDIMENTS            ! 
      !      (WITHOUT CONSOLIDATION: UNIFORM SEDIMENT BED) !                                   ! 
      ! ******************************************* *****  !

      IF(.NOT.TASS) THEN 
                 
        DO I = 1, NPOIN
          USTARP =SQRT(TAUP%R(I)/XMVE)
          IF(VITCE.GT.1.D-8) THEN
            AUX = MAX(((USTARP/VITCE)**2 - 1.D0),0.D0)
          ELSE
            AUX = 0.D0
          ENDIF
          FLUER%R(I) = PARTHENIADES*AUX
        ENDDO
              
      ELSE
 
      ! **************************************************** !
      ! Ib - FORMULATION FOR COHESIVE SEDIMENTS  + TASSEMENT ! 
      !      (WITH BEDLOAD)                                  ! 
      ! **************************************************** !
      
        DO I=1,NPOIN
C           
          DO J=1,NCOUCH_TASS
            IF(TAUP%R(I).GT.TOCE_VASE(J))THEN
              FLUER_LOC(J)=PARTHENIADES*
     &              ((TAUP%R(I)/TOCE_VASE(J))-1.D0)
            ELSE
              FLUER_LOC(J)=0.D0
            ENDIF
          ENDDO
          QER_VASE = 0.D0 
          TEMPS= DT
C          
          DO J= 1, NCOUCH_TASS 
            IF(ES(I,J).GE.1.D-6) THEN
C             CALCUL DE LA MASSE POTENTIELLEMENT ERODABLE DANS LA COUCHE J (KG/M2)      
              QE_COUCHE = FLUER_LOC(J) *XMVS * TEMPS          
              IF(QE_COUCHE.LT.MS_VASE(I,J)) THEN
                QER_VASE = QER_VASE  + QE_COUCHE
                GO TO 10             
              ELSE
                QER_VASE = QER_VASE + MS_VASE(I,J)
                TEMPS= TEMPS-MS_VASE(I,J)/FLUER_LOC(J)/XMVS
                TEMPS=MAX(TEMPS,0.D0)
              ENDIF
            ENDIF
          ENDDO
C
          IF(LNG.EQ.1) THEN
            WRITE(LU,*) 'ATTENTION TOUTES LES COUCHES SONT VIDES'
          ENDIF
          IF(LNG.EQ.2) THEN
            WRITE(LU,*) 'BEWARE, ALL LAYERS EMPTY'
          ENDIF
          CALL PLANTE(1)
          STOP           

10        CONTINUE   
C         attention partheniades est déja divise par XMVS?
          FLUER%R(I) = QER_VASE/DT/XMVS 
C
        ENDDO
      ENDIF
!
!======================================================================!
!======================================================================!
!
      RETURN      
      END
