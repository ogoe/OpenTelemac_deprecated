C                *************************
                 SUBROUTINE CONDIS_SISYPHE
C                *************************
C
     *(CONSTFLOW)
C
C**********************************************************************
C SISYPHE VERSION 5.5                             B. MINH DUC
C                                                 F. HUVELIN
C
C
C BUI MINH DUC AOUT 2003
C Creation de la subroutine 
C
C F. HUVELIN FEV. 2004 
C Ajout de la subroutine a sisyphe v5p4
C
C COPYRIGHT EDF-BAW
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
C |                |    |                                              C
C |                |    |                                              C
C |________________|____|______________________________________________C
C                                                                      C
C ---------------------------------------------------------------------C
C
C ORGANISATION
C ------------
!
! PROGRAMME APPELANT :
! --------------------
! SISYPHE
!
! PROGRAMMES APPELES : 
! --------------------
!
!=======================================================================!
!=======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS                !
!=======================================================================!
!=======================================================================!
!
!1/ MODULES
!----------
!
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_SISYPHE
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      LOGICAL, INTENT(INOUT) :: CONSTFLOW
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, EXTERNAL      :: P_ISUM
!
!3/ VARIABLES LOCALES
!--------------------
!
      INTEGER          :: NZFMAX, I
      DOUBLE PRECISION :: ZFMAX, C
!
!=======================================================================!
!=======================================================================!
!                               PROGRAMME                               !
!=======================================================================!
!=======================================================================!
!
      NZFMAX = 0
!
      IF(CONSTFLOW) THEN
         CALL OS('X=X+Y   ', ECPL, E, S, C)
!  
         DO I=1,NPOIN
            ZFMAX = ABS(ECPL%R(I)) - CRIT_CFD*HCPL%R(I)
            IF (ZFMAX.GT.1.D-8) NZFMAX=NZFMAX+1
         ENDDO
!  
         IF (NCSIZE.GT.1) THEN
            NZFMAX=P_ISUM(NZFMAX)
            CALL PARCOM(ECPL,2,MESH)
         ENDIF
!
         IF (NZFMAX.GE.1) CONSTFLOW = .FALSE.
      ENDIF
!
      IF(.NOT.CONSTFLOW) THEN
         CALL OS('X=C     ', ECPL,  S, S, 0.D0)
         CALL OS('X=Y     ', HCPL, HN, S,    C)
!
         IF (NCSIZE.GT.1) THEN
            CALL PARCOM(ECPL,2,MESH)
            CALL PARCOM(HCPL,2,MESH)
         ENDIF
      ENDIF
!
!=======================================================================!
!=======================================================================!
!
      RETURN
      END
