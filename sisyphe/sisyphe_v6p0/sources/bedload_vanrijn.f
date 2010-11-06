      ! ************************** !
        SUBROUTINE BEDLOAD_VANRIJN ! (_IMP_)
      ! ************************** !

     &  (TOB,MU, NPOIN, DM, DENS, GRAV, DSTAR, AC, QSC)


C**********************************************************************C
C SISYPHE VERSION 5.4       2004  C. VILLARET                          C
C SISYPHE VERSION 5.2  Oct. 2001  BUI MINH DUC                         C
C**********************************************************************C


               ! ====================================== !
               ! Bed-load transport formula of Van Rijn !
               ! ====================================== !


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
! CALL      ------                                                     !
!                                                                      !
!======================================================================!
!======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS               !
!======================================================================!
!======================================================================!

      ! 1/ MODULES
      ! ----------
      USE INTERFACE_SISYPHE,EX_BEDLOAD_VANRIJN => BEDLOAD_VANRIJN
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU                                                                          

      ! 2/ GLOBAL VARIABLES
      ! -------------------
      TYPE(BIEF_OBJ),   INTENT(IN)  :: TOB,MU
      INTEGER,          INTENT(IN)  :: NPOIN
      DOUBLE PRECISION, INTENT(IN)  :: DM, DENS, GRAV, DSTAR, AC
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: QSC


      ! 3/ LOCAL VARIABLES
      ! ------------------
      INTEGER          :: I
      DOUBLE PRECISION :: C1, C2, TETAP, T


!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!

      C1 = DENS * GRAV * DM
      C2 = 0.053D0 * SQRT(DM**3.0D0*DENS*GRAV) * DSTAR**(-0.3D0)

      DO I = 1, NPOIN

         ! ********************** !
         ! I -ADIM. SKIN FRICTION ! (_IMP_)
         ! ********************** !
         TETAP = 0.5D0 * TOB%R(I) * MU%R(I)**2 / C1


         ! ****************************** !
         ! II - TRANSPORT STAGE PARAMETER ! (_IMP_)
         ! ****************************** !
         IF (TETAP <= AC) THEN
            T=0.D0
         ELSE           
            T = (TETAP-AC)/MAX(AC,1.D-06)
         ENDIF


         ! ***************************** !
         ! III - BED LOAD TRANSPORT RATE ! (_IMP_)
         ! ***************************** !
         QSC%R(I) = C2 * (T**2.1D0)

       ENDDO

!======================================================================!
!======================================================================!

      RETURN
      END SUBROUTINE BEDLOAD_VANRIJN
