C                ******************
                 SUBROUTINE GF_USER
C                ******************
C
     *(TBEG_GF,TEND_GF,AT0)
C
C**********************************************************************
C SISYPHE VERSION 5.7                             B. MINH DUC
C                                                 F. HUVELIN

      ! ========================================== !
      ! FOR SUPPLIED BED-MATERIALS (grain feeding) !
      ! ========================================== !

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
C | NPOIN          | => |Mesh nodes number                             C
C | NSICLA         | => |Number of size-classes of bed material        C
C | X              | => |X coordinate of the mesh                      C
C | Y              | => |Y coordinate of the mesh                      C
C | AT0            | => |Time of the simulation                        C
C | FRACSED        | <=>|Fraction of size-classes of bed material      C
C                       (default value = initial value of case file)   C
C | TBEG           | <= |Time of feeding-begin                         C
C | TEND           | <= |Time of feeding-end                           C
C | DZF_GF         | <= |Bed-level change                              C
C |________________|____|______________________________________________C
C                    <=  Can't be change by the user                   C
C                    =>  Can be changed by the user                    C 
C ---------------------------------------------------------------------C

!=======================================================================!
!=======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS                !
!=======================================================================!
!=======================================================================!
!
      USE DECLARATIONS_SISYPHE
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION, INTENT(INOUT) :: TBEG_GF, TEND_GF
      DOUBLE PRECISION, INTENT(IN)    :: AT0
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C 
      !2/ VARIABLES LOCALES
      !--------------------
!     INTEGER :: I
!     DOUBLE PRECISION :: X1,X2,X3,X4
!     DOUBLE PRECISION :: Y1,Y2,Y3,Y4
!     DOUBLE PRECISION, POINTER :: X(:),Y(:)
!     X=>MESH%X%R
!     Y=>MESH%Y%R

!=======================================================================!
!=======================================================================!
!                               PROGRAMME                               !
!=======================================================================!
!=======================================================================!

! EXAMPLE : (to be filled by the user)
! ---------
C      TBEG_GF = 1.D0
C      TEND_GF = 2.D0
C
C      X1 = 0.D0
C      Y1 = 0.D0
C
C      X2 = 1.D0
C      Y2 = 0.D0
C
C      X3 = 0.D0
C      Y3 = 1.D0
C
C      X4 = 1.D0
C      Y4 = 1.D0
C
C      IF (AT0.GE.TBEG_GF.AND.AT0.LE.TEND_GF) THEN
C         DO I = 1, NPOIN
C            IF ((X(I)>X1).AND.(Y(I)>Y1).AND.
C     &          (X(I)<X2).AND.(Y(I)>Y2).AND.
C     &          (X(I)<X3).AND.(Y(I)<Y3).AND.
C     &          (X(I)>X4).AND.(Y(I)<Y4)     ) THEN
C                DZF_GF%R(I) = 1.25D-5
C            ENDIF
C         ENDDO
C      ENDIF
C
C
C  FOLLOWING LINES NEED TO BE COMMENTED OUT
C
      IF(LNG.EQ.1) WRITE(LU,52)
      IF(LNG.EQ.2) WRITE(LU,53)
C
52    FORMAT(/,1X,' STOP :',/
     *     ,1X,' GF_USER DOIT ETRE PROGRAMMEE')
53    FORMAT(/,1X,'SISYPHE IS STOPPED : ',/
     *      ,1X,' GF_USER MUSt BE IMPLEMENTED')
      CALL PLANTE(1)
      STOP
!
!=======================================================================
!
      RETURN
      END
