      ! ***************************** !
        SUBROUTINE SUSPENSION_FREDSOE ! 
      ! ***************************** !

     &  (ACLADM, TAUP, NPOIN, GRAV, 
     &   XMVE, XMVS, ZERO, AC,  CSTAEQ)


C**********************************************************************C
C SISYPHE VERSION 5.6  04/01/05  F. HUVELIN                            C
C SISYPHE VERSION 5.5  14/04/04  C. VILLARET  01 30 87 83 28           C
C**********************************************************************C


         ! ==================================================== !
         !   Reference concentration calculation at z= 2*d50    !
         ! thanks to the formula of Zyserman and Fredsoe (1994) !
         ! ==================================================== !

C
C
C 13/06/2008 : JMH : OPTIMISATION FORMULE AVEC AUX
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
C |   ACLADM       | => |
C |   CF           | => |
C |   TOB          | => |
C |   HCLIP        | => |
C |   AVA          | => |
C |   NPOIN        | => |
C |   CHARR        | => |
C |   KSPRATIO     | => |
C |   HMIN         | => |
C |   GRAV         | => |
C |   XMVE         | => |
C |   XMVS         | => |
C |   AC           | <=>|
C |   FLUER        | <= |
C !________________|____|______________________________________________C
C                    <=  Can't be change by the user                   C
C                    =>  Can be changed by the user                    C 
C ---------------------------------------------------------------------C
!                                                                      !
! CALLED BY SUSPENSION_FLUX                                            !
!                                                                      !
! CALL      ------                                                     !
!                                                                      !
!======================================================================!
!======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS               !
!======================================================================!
!======================================================================!

      ! 1/ MODULES
      ! ----------
      USE INTERFACE_SISYPHE,EX_SUSPENSION_FREDSOE => SUSPENSION_FREDSOE
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU


      ! 2/ GLOBAL VARIABLES
      ! -------------------
      TYPE(BIEF_OBJ),   INTENT(IN)    :: ACLADM, TAUP
      INTEGER,          INTENT(IN)    :: NPOIN
      DOUBLE PRECISION, INTENT(IN)    :: GRAV, XMVE, XMVS
      DOUBLE PRECISION, INTENT(IN)    :: ZERO,AC
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: CSTAEQ


      ! 3/ LOCAL VARIABLES
      ! ------------------
      INTEGER            :: I
      DOUBLE PRECISION   ::  TETAP,AUX
C 
      DOUBLE PRECISION   :: CMAX
C
C     MAXIMUM CONCENTRATION CORRESPONDING TO DENSE PACKING
C
      DATA CMAX/0.6D0/
      INTRINSIC MAX

!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!

      ! ******************************** !
      !    I - CRITICAL SHIELD PARAMETER ! 
      ! ******************************** !    

      DO I=1,NPOIN

         ! ****************** !
         ! II - SKIN FRICTION ! 
         ! ****************** !                
         
         TETAP = TAUP%R(I) / (GRAV*(XMVS-XMVE)*ACLADM%R(I))

         ! ***************** !
         ! IV - EROSION FLUX ! (_IMP_)
         ! ***************** !
         ! Concentration increased by AVA because it is assumed 
         ! that it is computed only with one class of sediment
         
         IF(TETAP.GT.AC) THEN
           AUX=(TETAP-AC)**1.75D0
           CSTAEQ%R(I) = 0.331D0*AUX/(1.D0+0.72D0*AUX)
           CSTAEQ%R(I) = MIN(CSTAEQ%R(I),CMAX)
         ELSE
           CSTAEQ%R(I) = 0.D0
         ENDIF

      ENDDO

!======================================================================!
!======================================================================!

      RETURN      
      END SUBROUTINE SUSPENSION_FREDSOE
