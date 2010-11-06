      ! ************************** !
        SUBROUTINE BEDLOAD_BAILARD !
      ! ************************** !

     &(U2D,V2D,UCMOY,TOB,TOBW,THETAW,UW,FW,CF,NPOIN,PI,
     & XMVE,GRAV,DENS,XWC,ALPHAW,QSCX,QSCY,QSSX,QSSY,
     & UC3X,UC3Y,US4X,US4Y,THETAC,FCW,QSC,QSS,HOULE)


C**********************************************************************
C SISYPHE VERSION 6.0  01/10/2003                    C. VILLARET (LNHE)
C
C JMH LE 21/12/2006: BEDLOAD_CALCBAIL SUPPRIME, ARGUMENT HOULE ADDED                 
C**********************************************************************


                       ! ================== !
                       ! Formula of Bailard !
                       ! ================== !

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
! CALLED BY BEDLOAD_FORMULA                                            !
!                                                                      !
! CALL      BEDLOAD_DIRECTION                                          !
!           BEDLOAD_CALCBAIL                                           !
!           BEDLOAD_INTERACT                                           !
!                                                                      !
!======================================================================!
!======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS               !
!======================================================================!
!======================================================================!
!
! 1/ MODULES
! 
      USE INTERFACE_SISYPHE,EX_BEDLOAD_BAILARD => BEDLOAD_BAILARD
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
! 2/ GLOBAL VARIABLES
! 
      TYPE(BIEF_OBJ),   INTENT(IN)    :: U2D,V2D,UCMOY, TOB
      TYPE(BIEF_OBJ),   INTENT(IN)    :: TOBW, THETAW, UW, FW, CF
      INTEGER,          INTENT(IN)    :: NPOIN
      LOGICAL,          INTENT(IN)    :: HOULE
      DOUBLE PRECISION, INTENT(IN)    :: PI, XMVE, GRAV, DENS, XWC
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: ALPHAW        ! work array T1 
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: QSCX, QSCY    ! work array T2 and T3
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: QSSX, QSSY    ! work array T4 and T5
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: UC3X, UC3Y    ! work array T6 and T7
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: US4X, US4Y    ! work array T8 and T9
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: THETAC, FCW   ! work array T10 and T11
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: QSC, QSS
!
! 3/ LOCAL VARIABLES
! 
      INTEGER                     :: I
      DOUBLE PRECISION            :: C3, C4, PHI
      DOUBLE PRECISION, PARAMETER :: EPSC = 0.21D0   ! efficac. du char.
      DOUBLE PRECISION, PARAMETER :: EPSS = 0.025D0  ! efficac. de la susp.
      DOUBLE PRECISION            :: U3X, U3Y, NUM
!
!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!
!
!     CASE WITH WAVES
!
      IF(HOULE) THEN
!
!     ANGLE OF VELOCITY WITH Ox (IN RADIANS)
!
      CALL BEDLOAD_DIRECTION(U2D,V2D,NPOIN,PI,THETAC)
!
!     ANGLE OF WAVES WITH Ox (IN RADIANS)
!
      CALL OS('X=CY    ', X=ALPHAW, Y=THETAW, C=-PI/180.D0)
      CALL OS('X=X+C   ', X=ALPHAW, C=0.5D0*PI)
      CALL OS('X=Y-Z   ', X=ALPHAW, Y=ALPHAW, Z=THETAC)
!
!     PARAMETERS <|U|^2.U> , <|U|^3>  
!
!    
!     US4X AND US4Y ARE WORK ARRAYS, THEIR STRUCTURE IS GIVEN HERE
!     THE STRUCTURE OF THETAC (CATHERINE DON'T REMOVE THIS PLEASE)
      CALL CPSTVC(THETAC,US4X)
      CALL CPSTVC(THETAC,US4Y)
!      
      DO I = 1, NPOIN

         ! ********************* !
         ! I - REPERE DU COURANT !
         ! ********************* !
         
         U3X = UCMOY%R(I)**3
     &       + UCMOY%R(I)*UW%R(I)**2 * (1 + COS(2.D0*ALPHAW%R(I))/ 2.D0)
         U3Y = UCMOY%R(I)*UW%R(I)**2 * SIN(2.D0*ALPHAW%R(I)) / 2.D0

         ! ********************************************** !
         ! II - MOMENTS D'ORDRE 3 (HYP. DE HOUE LINEAIRE) ! 
         ! ********************************************** !
         
         UC3X%R(I) = U3X * COS(THETAC%R(I)) - U3Y * SIN(THETAC%R(I))
         UC3Y%R(I) = U3X * SIN(THETAC%R(I)) + U3Y * COS(THETAC%R(I))      

         ! ************************************************************ !
         ! III -  MOMENTS D'ORDRE 4 (HYP. DE HOULE ET COURANTS ALIGNES) ! 
         ! ************************************************************ !
         
         NUM = ( 8.D0*UCMOY%R(I)**4 + 3.D0*UW%R(I)**4
     &           + 24.D0*(UCMOY%R(I)**2)*(UW%R(I)**2) )*0.125D0
         US4X%R(I) = NUM * COS(THETAC%R(I))
         US4Y%R(I) = NUM * SIN(THETAC%R(I))

       ENDDO

      ! *********************************************** !
      ! IV -  COEFFICIENT DE FROTTEMENT HOULE + COURANT !
      ! *********************************************** !
      
      CALL BEDLOAD_INTERACT
     &     (UCMOY,TOBW,TOB,ALPHAW,FW,CF,UW,NPOIN,XMVE,FCW)

      ! ******************************** !
      ! V - CALCUL DES TAUX DE TRANSPORT !
      ! ******************************** !
      
      PHI = PI   / 6.D0  ! angle de frottement
      C3  = EPSC / (GRAV*DENS*TAN(PHI))
      C4  = EPSS / (GRAV*DENS*XWC)
      CALL OS('X=CYZ   ', X=QSCX, Y=FCW,  Z=UC3X, C=C3)
      CALL OS('X=CYZ   ', X=QSCY, Y=FCW,  Z=UC3Y, C=C3)
      CALL OS('X=CYZ   ', X=QSSX, Y=FCW,  Z=US4X, C=C4)
      CALL OS('X=CYZ   ', X=QSSY, Y=FCW,  Z=US4Y, C=C4)
!
!     CASE WITHOUT WAVES
!
      ELSE
!
        CALL PLANTE(1)
        STOP 'BAILARD WITHOUT WAVES NOT PROGRAMMED'
!
      ENDIF
!
!     NORMS
!
      CALL OS('X=N(Y,Z)', X=QSC,  Y=QSCX, Z=QSCY)
      CALL OS('X=N(Y,Z)', X=QSS,  Y=QSSX, Z=QSSY)

!======================================================================!
!======================================================================!

      RETURN
      END
