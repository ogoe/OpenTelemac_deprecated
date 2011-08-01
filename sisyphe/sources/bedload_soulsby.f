      ! ************************** !
        SUBROUTINE BEDLOAD_SOULSBY ! (_IMP_)
      ! ************************** !

     &  (UCMOY,HN, UW, NPOIN, DENS, GRAV, DM, DSTAR, HMIN, D90, QSC,
     &   QSS)

C**********************************************************************C
C SISYPHE VERSION 5.4  Nov. 2003   C.VILLARET                          C
C SISYPHE VERSION 5.2  22/05/2001  SOGREAH                             C
C**********************************************************************C


           ! ================================================ !
           ! Bed-load transport formula of Soulsby & Van Rijn !
           ! ================================================ !


C COPYRIGHT EDF-BAW-IFH
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
C |________________|____|______________________________________________C
C                    <=  Can't be change by the user                   C
C                    =>  Can be changed by the user                    C 
C ---------------------------------------------------------------------C
!                                                                      !
! CALLED BY BEDLOAD_SOLIDISCHARGE                                      !
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
      USE INTERFACE_SISYPHE,EX_BEDLOAD_SOULSBY => BEDLOAD_SOULSBY
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU

      ! 2/ GLOBAL VARIABLES
      ! -------------------
      TYPE(BIEF_OBJ),   INTENT(IN)  :: HN, UCMOY, UW
      INTEGER,          INTENT(IN)  :: NPOIN
      DOUBLE PRECISION, INTENT(IN)  :: DENS, GRAV, DM, DSTAR, HMIN, D90
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: QSC, QSS


      ! 3/ LOCAL VARIABLES
      ! ------------------
      INTEGER                     :: I
      DOUBLE PRECISION            :: COEF, ASS, ASB, CD
      DOUBLE PRECISION            :: UCR, VTOT, TRA
      DOUBLE PRECISION, PARAMETER :: Z0=0.006D0


!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!

      ! ************************* !
      ! I - COEFFCIENT SUSPENSION ! 
      ! ************************* !
      COEF = (DENS *GRAV*DM)**1.2D0
      ASS  = 0.012D0*DM*(DSTAR**(-0.6D0))/COEF


      DO I = 1, NPOIN
        
         
         ! *************************** !
         ! III - COEFFICIENT CHARRIAGE ! 
         ! *************************** !
         ASB = 0.005D0*HN%R(I)*(DM/MAX(HN%R(I),DM))**1.2D0 / COEF


         ! ********************************** !
         ! IV - COEFFICIENT DE RUGOSITE CD    ! 
         !      soulsby: z0=0.006 --> ks=18cm ! 
         ! ********************************** !
         CD = (0.4D0 / (LOG(MAX(HN%R(I),Z0)/Z0)-1.D0))**2


         ! ************************************************ !
         ! V - CALCUL DE LA VITESSE DE COURANT CRITIQUE UCR ! 
         ! ************************************************ !
         IF (DM < 0.0005D0) THEN
            UCR = 0.19D0*(DM**0.1D0)*LOG10(4.D0*MAX(HN%R(I),D90)/D90)
         ELSE
            UCR = 8.50D0*(DM**0.6D0)*LOG10(4.D0*MAX(HN%R(I),D90)/D90)
         ENDIF


         ! ************************************************* !
         ! VI - VITESSE INDUITE PAR LE COURANT ET LES VAGUES ! 
         ! ************************************************* !
         VTOT = SQRT(UCMOY%R(I)**2+(0.018D0/CD)*UW%R(I)**2)


         ! *********************************************** !
         ! VII - TRANSPORT PAR SUSPENSION ET PAR CHARRIAGE ! 
         ! *********************************************** !
         IF (VTOT > UCR) THEN
            TRA     = UCMOY%R(I)  * (VTOT - UCR )**2.4D0
            QSS%R(I)= ASS * TRA
            QSC%R(I)= ASB * TRA
         ELSE
            QSS%R(I) = 0.D0
            QSC%R(I) = 0.D0
         ENDIF
      ENDDO

!======================================================================!
!======================================================================!

      RETURN
      END SUBROUTINE BEDLOAD_SOULSBY
