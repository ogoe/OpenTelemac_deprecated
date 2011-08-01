C                       ************************
                        SUBROUTINE FRICTION_INIT
C                       ************************
C
C***********************************************************************
C  TELEMAC-2D VERSION 6.0             J-M HERVOUET (LNHE) 01 30 87 80 18
C***********************************************************************
C
C 20/04/04 : subroutine written by F. Huvelin
C
C
C         ! ----------------------------------------------- !
C         !   Friction calculation by zone initialization   !
C         ! ----------------------------------------------- !
C
C
C               TTTTT EEEEE L     EEEEE M   M   AA  CCCCC
C                 T   E     L     E     MM MM  A  A C
C                 T   EEE   L     EEE   M M M  AAAA C
C                 T   E     L     E     M   M  A  A C
C                 T   EEEEE LLLLL EEEEE M   M  A  A CCCCC
C
C
C----------------------------------------------------------------------C
C                             ARGUMENTS                                C
C .________________.____.______________________________________________C
C |      NOM       |MODE|                   ROLE                       C
C |________________|____|______________________________________________C
C |                | => |                                              C
C |________________|____|______________________________________________C
C                    <=  input value                                   C
C                    =>  output value                                  C 
C ---------------------------------------------------------------------C
!
!=======================================================================!
!=======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS                !
!=======================================================================!
!=======================================================================!
!
      USE BIEF
      USE FRICTION_DEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_TELEMAC2D
      USE INTERFACE_TELEMAC2D
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
      !2/ Local variables
      !------------------
      INTEGER :: I, J, K
      LOGICAL :: FRICTION_ERR
!
!=======================================================================!
!=======================================================================!
!                               PROGRAMME                               !
!=======================================================================!
!=======================================================================!
!
      CALL FRICTION_READ(T2D_FILES(T2DCOF)%LU,
     *                   NZONMX,ITURB,LISRUG,LINDNER,
     *                   T2D_FILES(T2DCOF)%NAME,NZONES,FRTAB)
!
      ! Initialization (all elements with -1
      ! in order to check after user initialization)
      ! --------------------------------------------
      DO I = 1, CF%DIM1
        KFROPT%I(I) = -1
      ENDDO
!
      ! User initialization
      ! -------------------
      CALL FRICTION_USER
!         
      FRICTION_ERR = .FALSE.
!
      ! Check value
      ! -----------
! FH : FOR QUASI-BUBBLE
! FH : 2004/03/01
! =>
      DO I=1, NPOIN
! <=
! FH : 2004/03/01
! FH : FOR QUASI-BUBBLE
            
         ! No friction zone defined
         ! ------------------------
         IF (KFROPT%I(I) == -1) THEN

            FRICTION_ERR = .TRUE.
            IF(NCSIZE>1) THEN
              K = MESH%KNOLG%I(I)
            ELSE
              K = I
            ENDIF
            IF(LNG == 1) WRITE(LU,10) K
            IF(LNG == 2) WRITE(LU,11) K
!
         ! Local numbering of the zone
         ! ---------------------------
         ELSE
            DO J = 1, NZONES
               IF(KFROPT%I(I) == FRTAB%ADR(J)%P%GNUMB(1)) THEN
                  KFROPT%I(I) = J
                  EXIT
               ENDIF
               IF(J==NZONES) THEN
                  FRICTION_ERR = .TRUE.
                  IF (NCSIZE>1) THEN
                      K = MESH%KNOLG%I(I)
                  ELSE
                      K=I
                  ENDIF
                  IF (LNG==1) WRITE(LU,20) K,KFROPT%I(I)
                  IF (LNG==2) WRITE(LU,21) K,KFROPT%I(I)
               ENDIF
            ENDDO
         ENDIF
!
      ENDDO
!
10    FORMAT('AUCUNE ZONE DE FROTTEMENT DEFINI POUR LE NOEUD : ',i5)
11    FORMAT('NO FRICTION ZONE DEFINED FOR THE NODE : ',i5)
!
20    FORMAT('MAUVAISE INITIALISATION DE LA ZONE DE FROTTEMENT DU ',
     &       'NOEUD :',i5
     &     ,/' ZONE : ',i9,' INCONNUE')
21    FORMAT('WRONG INITIALIZATION OF THE FRICTION ZONE FOR THE NODE :'
     &     ,  i5
     &     ,/' ZONE : ',i9,' UNKNOWN')
!
      IF(FRICTION_ERR) THEN
        CALL PLANTE(1)
        STOP
      ENDIF
!
! FH : FOR QUASI-BUBBLE
! FH : 2004/03/01
! =>
      ! Vector initialization : whole domain
      ! (for quasi_bubble, see friction_choice.f : CALL FRICTION_BUBBLE)
      ! ----------------------------------------------------------------
      IF (LINDNER) THEN
         DO I = 1, NPOIN
            CHESTR%R(I) = FRTAB%ADR(KFROPT%I(I))%P%RCOEF(1)
            NDEFMA%R(I) = FRTAB%ADR(KFROPT%I(I))%P%NDEF (1)
            NKFROT%I(I) = FRTAB%ADR(KFROPT%I(I))%P%RTYPE(1)
            LINDDP%R(I) = FRTAB%ADR(KFROPT%I(I))%P%DP
            LINDSP%R(I) = FRTAB%ADR(KFROPT%I(I))%P%SP
         ENDDO
      ELSE
         DO I = 1, NPOIN
            CHESTR%R(I) = FRTAB%ADR(KFROPT%I(I))%P%RCOEF(1)
            NDEFMA%R(I) = FRTAB%ADR(KFROPT%I(I))%P%NDEF (1)
            NKFROT%I(I) = FRTAB%ADR(KFROPT%I(I))%P%RTYPE(1)
         ENDDO
      ENDIF
!
      ! Vector initialization : boundary conditions
      ! -------------------------------------------
      IF ((ITURB == 3).AND.(LISRUG == 2)) THEN
         DO J = 1, MESH%NPTFR
            I = MESH%NBOR%I(J)
            CHBORD%R(J) = FRTAB%ADR(KFROPT%I(I))%P%RCOEF(2)
            NDEF_B%R(J) = FRTAB%ADR(KFROPT%I(I))%P%NDEF (2)
            KFRO_B%I(J) = FRTAB%ADR(KFROPT%I(I))%P%RTYPE(2)
         ENDDO
      ENDIF
!
! <=
! FH : 2004/03/01
! FH : FOR QUASI-BUBBLE
      ! KFROT is used in order to know
      ! how many zone haves a friction coeffcient
      ! -----------------------------------------
      KFROT = 0
      DO I =1, NZONES
        IF(FRTAB%ADR(I)%P%RTYPE(1).NE.0) KFROT = KFROT + 1
      ENDDO
!
!=======================================================================!
!=======================================================================!
!
      RETURN
      END SUBROUTINE FRICTION_INIT
