!       ********************************
        SUBROUTINE BEDLOAD_HIDING_FACTOR
!       ********************************
!
     &(ACLADM, HIDFAC, NPOIN, HIDI, DM, KARIM_HOLLY_YANG, HIDING)
!
!**********************************************************************C
! SISYPHE VERSION 5.5  14/09/2004  F. HUVELIN                          C
! SISYPHE VERSION 5.3  --/--/2002  M. GONZALES DE LINARES              C
! SISYPHE VERSION 5.3  --/11/2002  B. MINH DUC                         C
!**********************************************************************C
!
!
!     ! ========================================================= !
!     ! Hiding factor for each node, sediment class and time step !
!     ! ========================================================= !
!
!
! COPYRIGHT EDF-BAW-IFH   
!**********************************************************************C
!                                                                      C
!                 SSSS I   SSSS Y   Y PPPP  H   H EEEEE                C
!                S     I  S      Y Y  P   P H   H E                    C
!                 SSS  I   SSS    Y   PPPP  HHHHH EEEE                 C
!                    S I      S   Y   P     H   H E                    C
!                SSSS  I  SSSS    Y   P     H   H EEEEE                C
!                                                                      C
!----------------------------------------------------------------------C
!                             ARGUMENTS                                C
! .________________.____.______________________________________________C
! |      NOM       |MODE|                   ROLE                       C
! |________________|____|______________________________________________C
! |________________|____|______________________________________________C
!
! ---------------------------------------------------------------------C
!                                                                      !
! CALLED BY SISYPHE                                                    !
!                                                                      !
! CALL                                                                 !
!                                                                      !
!======================================================================!
!======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS               !
!======================================================================!
!======================================================================!
!
! 1/ MODULES
! ----------
!
      USE INTERFACE_SISYPHE,
     &    EX_BEDLOAD_HIDING_FACTOR => BEDLOAD_HIDING_FACTOR
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!      
!
! 2/ GLOBAL VARIABLES
! -------------------
!
      TYPE(BIEF_OBJ),   INTENT(IN)  :: ACLADM
      INTEGER,          INTENT(IN)  :: HIDFAC, NPOIN
      DOUBLE PRECISION, INTENT(IN)  :: HIDI, DM, KARIM_HOLLY_YANG
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: HIDING
!
!
! 3/ LOCAL VARIABLES
! ------------------
!
      INTEGER          :: J
      DOUBLE PRECISION :: C1, C2
!
!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!
!
! *************************** !
! Ia - CONSTANT HIDING FACTOR !
! *************************** !
!
      IF (HIDFAC == 0) THEN
!
         CALL OS('X=C     ', X=HIDING, C=HIDI)
!
! ************************** !
! Ib - FORMULA OF EGIAZAROFF ! 
! ************************** !
!
      ELSEIF (HIDFAC == 1) THEN
!
         C1 = LOG10(19.D0)
         C2 = 19.D0*DM
         DO J = 1, NPOIN
           HIDING%R(J) = (C1/LOG10(C2/ACLADM%R(J)))**2
         ENDDO
!
! ********************************** !
! Ic - FORMULA OF ASHIDA AND MICHIUE !
! ********************************** !
!
      ELSEIF (HIDFAC == 2) THEN
!
         C1 = LOG10(19.D0)
         C2 = 19.D0*DM
         DO J = 1, NPOIN
!
            IF(DM/ACLADM%R(J) >= 0.4D0) THEN
              HIDING%R(J) = (C1 / LOG10(C2/ACLADM%R(J)) )**2
            ELSE
              HIDING%R(J) = 0.85D0*(ACLADM%R(J)/DM)
            ENDIF
!
         ENDDO
!
! ************************************* !
! Ie - FORMULA OF KARIM, HOLLY AND YANG !
! ************************************* !
!      
      ELSEIF (HIDFAC == 4) THEN
!
         CALL OS('X=1/Y   ', X=HIDING, Y=ACLADM)
         CALL OS('X=CX    ', X=HIDING, C=DM)
         CALL OS('X=Y**C  ', X=HIDING, Y=HIDING, C=KARIM_HOLLY_YANG)
!
      ELSE
!
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'FORMULE DE MASQUAGE INCONNUE : ',HIDFAC
        ENDIF 
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) 'UNKNOWN HIDING FACTOR FORMULA: ',HIDFAC
        ENDIF 
        CALL PLANTE(1)
        STOP     
!      
      ENDIF
!
!======================================================================!
!======================================================================!
!
      RETURN
      END SUBROUTINE BEDLOAD_HIDING_FACTOR
