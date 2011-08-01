      ! ************************ !
        SUBROUTINE DEBUG_SISYPHE ! (_IMP_)
      ! ************************ !
     * (NAME, ILOOP, NLOOP)

C**********************************************************************C
C SISYPHE VERSION 5.5                           F. HUVELIN    22/11/04 C


               ! =============== !
               ! Sisyphe debuger !
               ! =============== !


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
C |________________|____|______________________________________________C
C                    <=  Can't be change by the user                   C
C                    =>  Can be changed by the user                    C 
C ---------------------------------------------------------------------C
!                                                                      !
! CALLED BY                                                            !
!                                                                      !
!                                                                      !
!======================================================================!
!======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS               !
!======================================================================!
!======================================================================!
!
      ! 1/ MODULES
      ! ----------
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_SISYPHE
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
      ! 2/ GLOBAL VARIABLES
      ! -------------------
      INTEGER, PARAMETER :: SIZE =100
      CHARACTER(LEN=SIZE)  , INTENT(IN)           :: NAME
      INTEGER,   INTENT(IN), OPTIONAL :: ILOOP, NLOOP
!
      INTEGER :: LONGUEUR
!
!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!
!
      LONGUEUR = LEN(NAME)
!      
      WRITE(LU,*) NAME
      IF (PRESENT(ILOOP).AND.PRESENT(NLOOP)) THEN
         WRITE(LU,*) ILOOP, NLOOP
      ENDIF
!     JMH 03/02/05
!     MAY BE ADDED AT BAW BUT THIS IS NOT STANDARD FORTRAN
!     CALL FLUSH(LU)
!
100   FORMAT('LOOP : ',I8,'/',I8)
!
!======================================================================!
!======================================================================!
!
      RETURN
      END SUBROUTINE DEBUG_SISYPHE
