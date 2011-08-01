      ! ************************* !
        SUBROUTINE BEDLOAD_BIJKER !
      ! ************************* !
     &  (TOBW,TOB,MU,KSP,KSR,HN,NPOIN,DM,DENS,XMVE,GRAV,XWC, 
     &   KARMAN,ZERO,T4,T7,T8,T9,QSC,QSS,BIJK,HOULE)


C**********************************************************************C
C SISYPHE VERSION 5.6  Dec 2004  F. HUVELIN                            C
C SISYPHE VERSION 5.4  10/03/04  C. VILLARET                           C
C SISYPHE VERSION 5.1  26/11/01  E. BEN SLAMA                          C
C SISYPHE VERSION 5.1  26/11/01  T. BOULET                             C
C SISYPHE VERSION 5.1  26/11/01  C. MACHET                             C
C**********************************************************************C


              ! ==================================== !
              ! Bed-load transport formula of Bijker !
              ! ==================================== !


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
      USE INTERFACE_SISYPHE,EX_BEDLOAD_BIJKER => BEDLOAD_BIJKER
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU

      ! 2/ GLOBAL VARIABLES
      ! -------------------
      TYPE(BIEF_OBJ),   INTENT(IN)    :: TOBW, TOB, KSR,KSP, HN,MU
      INTEGER,          INTENT(IN)    :: NPOIN
      LOGICAL,          INTENT(IN)    :: HOULE
      DOUBLE PRECISION, INTENT(IN)    :: DM, DENS, XMVE, GRAV, XWC
      DOUBLE PRECISION, INTENT(IN)    :: KARMAN, ZERO
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: T4
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: T7, T8, T9
      TYPE(BIEF_OBJ),   INTENT(INOUT)   :: QSC, QSS


      ! 3/ LOCAL VARIABLES
      ! ------------------
      INTEGER                      :: I
      DOUBLE PRECISION             :: C1, C2, UCF
      DOUBLE PRECISION, INTENT(IN) :: BIJK

!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!

      ! ***************************************************** !
      ! I - CONTRAINTE SOUS L'ACTION COMBINE HOULE ET COURANT !
      ! ***************************************************** !
      
      IF(HOULE) THEN
        CALL OS('X=CY    ', X=T4, Y=TOBW, C= 0.5D0)
        CALL OS('X=X+Y   ', X=T4, Y=TOB)
      ELSE
        CALL OS('X=Y     ', X=T4, Y=TOB)
      ENDIF

      ! ******************************************************* !
      ! II - CORRECTION POUR PRISE EN COMPTE DES FORMES DE FOND !
      ! ******************************************************* !
      
C      CALL OS('X=Y/Z   ', X=MU, Y=CFP, Z=CF)
C      CALL OS('X=Y**C  ', X=MU, Y=MU , C=0.75D0)

      ! ***************************** !
      ! III - TRANSPORT PAR CHARRIAGE !
      ! ***************************** !
      C1 = BIJK*DM
      C2 = DENS*DM*XMVE*GRAV
      DO I = 1, NPOIN
         IF (T4%R(I)*MU%R(I)> ZERO) THEN
            QSC%R(I) = C1*SQRT(TOB%R(I)/XMVE )
     &               * EXP(-0.27D0*(C2/(T4%R(I)*MU%R(I))))
         ELSE
            QSC%R(I) = 0.D0
         ENDIF
      ENDDO

      ! *********************************************************** !
      ! IV- NOMBRE DE ROUSE ET BORNE INF. DE L'INTEGRALE D'EINSTEIN ! (_IMP_)
      ! *********************************************************** !
      DO I = 1, NPOIN
         IF (T4%R(I) > 0.D0) THEN                           
            UCF     = SQRT( T4%R(I) / XMVE) 
            T7%R(I) = XWC / ( KARMAN * UCF )
C            AUX     = 1.D0 + KARMAN*SQRT(2.D0/MAX(CF%R(I),ZERO))
C            T8%R(I) = 30.D0*EXP(-AUX)       
             T8%R(I) = MAX(KSR%R(I),KSP%R(I))/MAX(HN%R(I),ZERO)
         ELSE
            T7%R(I)= 100001.D0
            T8%R(I)= 100001.D0
         ENDIF
      ENDDO

      ! ************************************ !
      ! V - CALCUL DE L'INTEGRALE D'EINSTEIN ! (_IMP_)
      ! ************************************ !
      CALL INTEG(T7%R, T8%R, T9%R, NPOIN)

      ! ************************************** !
      ! VI - CALCUL DU TRANSPORT EN SUSPENSION ! (_IMP_)
      ! ************************************** !
      CALL OS('X=YZ    ', X=QSS, Y=T9, Z=QSC)

!======================================================================!
!======================================================================!

      RETURN
      END SUBROUTINE BEDLOAD_BIJKER
