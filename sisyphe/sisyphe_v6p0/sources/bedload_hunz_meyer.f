      ! ***************************** !
        SUBROUTINE BEDLOAD_HUNZ_MEYER ! 
      ! ***************************** !
C
     &  (TOB, MU, ACLADM, UNLADM, NPOIN, DENS, XMVE, GRAV, DM, AC,
     &    TETAP, AHUNZI, ACP, HIDING, QSC)
C
C**********************************************************************
C SISYPHE VERSION 6.0                        
C SISYPHE VERSION 5.4  Oct   2003  C. VILLARET                         
C SISYPHE VERSION 5.2  Jan.  2002  BUI MINH DUC                        
C**********************************************************************
C
C
C              ! ====================================== !
C              !   Bed-load formula of Hunziker (1995)  !
C              ! (adaptation of the Meyer-Peter formula !
C              ! ====================================== !
C
C
C   NOTE JMH : **-1.5 AND **1.5 SHOULD BE OPTIMIZED
C
C
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
C
C ---------------------------------------------------------------------C
!                                                                      !
! CALLED BY BEDLOAD_FORMULA                                            !
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
      USE INTERFACE_SISYPHE,
     &          EX_BEDLOAD_HUNZ_MEYER => BEDLOAD_HUNZ_MEYER
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU


      ! 2/ GLOBAL VARIABLES
      ! -------------------
      TYPE(BIEF_OBJ),   INTENT(IN)    :: TOB, MU, ACLADM, UNLADM
      INTEGER,          INTENT(IN)    :: NPOIN
      DOUBLE PRECISION, INTENT(IN)    :: DENS, XMVE, GRAV, DM, AC
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: TETAP, AHUNZI ! work array T1, T2
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: ACP           ! work array T3
      TYPE(BIEF_OBJ),   INTENT(INOUT)   :: HIDING, QSC


      ! 3/ LOCAL VARIABLES
      ! ------------------
      INTEGER          :: I
      DOUBLE PRECISION :: C1, C2

!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!

      ! ************************************* !
      ! I - CONTRAINTE DE PEAU ADIMENSIONELLE ! 
      ! ************************************* !
      C1 = 1.D0/(DENS*XMVE*GRAV*DM)
      CALL OS('X=CYZ   ', X=TETAP, Y=TOB, Z=MU, C=C1)
!
!     CHANGED BY JMH ON 28/10/2009 AFTER MODIFICATIONS BY
!     REBEKKA KOPMANN TRANSMITTED BY JACEK JANKOWSKI
!     CALL OS('X=+(Y,C)', X=TETAP , Y=TETAP, C= 1.D-06 )
      CALL OS('X=+(Y,C)', X=TETAP , Y=TETAP, C= 1.D-02 )
C
      CALL OS('X=Y**C  ', X=AHUNZI, Y=TETAP, C=-1.5D0  )
      CALL OS('X=CX    ', X=AHUNZI, C= 0.011D0)
      CALL OS('X=X+C   ', X=AHUNZI, C=-0.3D0  )
!
! RK COMMENT:
! AChtung AHUNZI kann so gross werden, dass der hiding factor unendlich wird
! Hunziker selbst schlaegt eine Begrenzung auf 2,3 vor.
! Wir haben erst mal eine Begrenzung bei ungefaehr 10 angenommen
! (ergibt sich aus tetap nicht kleiner als 0,01)
!
!     REMARK BY JMH: I WOULD STRONGLY RECOMMEND A SINGLE LOOP
!                    WITH THE WHOLE FORMULA, INSTEAD OF PILING
!                    UP CALLS TO OS
!
      DO I = 1, NPOIN
        HIDING%R(I) = (DM/ACLADM%R(I))**(-AHUNZI%R(I))
      ENDDO

      ! ************************************************* !
      ! IV - CORRECTION DE LA CONTRAINTE CRITIQUE ADIMENS ! 
      ! ************************************************* !
      CALL OS('X=Y/Z   ', X=ACP, Y=UNLADM, Z=ACLADM)
      CALL OS('X=Y**C  ', X=ACP, Y=ACP   , C=0.33D0)
      CALL OS('X=CX    ', X=ACP, C=AC)

      ! ********************* !
      ! V - TAUX DE TRANSPORT ! 
      ! ********************* !
      CALL OS('X=Y-Z   ', X=QSC, Y=TETAP , Z=ACP )
      CALL OS('X=+(Y,C)', X=QSC, Y=QSC   , C=0.D0)
      CALL OS('X=XY    ', X=QSC, Y=HIDING)
      CALL OS('X=Y**C  ', X=QSC, Y=QSC   , C=1.5D0)
      C2 = 5.D0*SQRT(GRAV*DENS*DM**3)
      CALL OS('X=CX    ', X=QSC, C=C2)

!======================================================================!
!======================================================================!

      RETURN
      END
