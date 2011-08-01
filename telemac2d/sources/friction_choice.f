C                       **************************
                        SUBROUTINE FRICTION_CHOICE
C                       **************************
C
     &(FRICTION_PASS,KARMAN)
C
C***********************************************************************
C  TELEMAC-2D VERSION 6.0             J-M HERVOUET (LNHE) 01 30 87 80 18
C***********************************************************************
C
C 20/04/04 : subroutine written by F. Huvelin
C
C
C         ! ----------------------------------------------- !
C         !   Main subroutine for the friction computation  !
C         ! ----------------------------------------------- !
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
C----------------------------------------------------------------------C
!
!======================================================================!
!======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS               !
!======================================================================!
!======================================================================!
!
      USE BIEF
      USE FRICTION_DEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_TELEMAC2D
      USE INTERFACE_TELEMAC2D, EX_FRICTION_CHOICE => FRICTION_CHOICE
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER,          INTENT(IN) :: FRICTION_PASS
      DOUBLE PRECISION, INTENT(IN) :: KARMAN
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER                     :: I
      DOUBLE PRECISION, PARAMETER :: VK = 1.D-6
!
!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!
!
! Initialization
! --------------
!
      IF(FRICTION_PASS == 0) THEN
!
! Zones initialization
! --------------------
!
         IF (FRICTB) THEN
            CALL FRICTION_INIT
            CALL STRCHE
! FH : FOR QUASI-BUBBLE
! FH : 2004/03/01
!jaj for quadratic elements 
! =>
            IF (CF%ELM /= H%ELM) THEN
              IF(CF%ELM==12 .AND. H%ELM==11) THEN
                CALL FRICTION_BUBBLE
     &              (IKLE, NPOIN, NELEM, NELMAX, LINDNER, NKFROT,
     &               CHESTR, NDEFMA, LINDDP, LINDSP)
              ELSEIF (CF%ELM==13 .AND. H%ELM==11) THEN
                CALL FRICTION_QUAD
     &              (IKLE%I, NPOIN, NELEM, NELMAX, LINDNER, NKFROT,
     &               CHESTR, NDEFMA, LINDDP, LINDSP)
!                WRITE (LU,*) 
!     &           'FRICTION_CHOICE::QUADRATIC ELEMENTS NOT IMPLEMENTED.'
!                CALL PLANTE(1)
              ELSE
                WRITE (LU,*)
     &           'FRICTION_CHOICE::DISCRETISATION NOT IMPLEMENTED.'
                WRITE (LU,*)
     &           'CF%ELM, H%ELM: ',CF%ELM, H%ELM
                CALL PLANTE(1)
              ENDIF
            ENDIF
! <=
!jaj for quadratic elements
! FH : 2004/03/01
! FH : FOR QUASI-BUBBLE

!
!        Uniform case
!        ------------
!
         ELSE

            ! Chestr for boundary conditions initialization
            ! -----------------------------------------------
            IF(KFROT.EQ.1.AND.ITURB.EQ.3.AND.LISRUG.EQ.2) THEN
              DO I = 1, NPTFR
                  CHBORD%R(I) = CHESTR%R(MESH%NBOR%I(I))
              ENDDO
            ELSEIF(ITURB.EQ.3.AND.LISRUG.EQ.2) THEN
              CALL OS('X=C     ', X=CHBORD, C=SB)
            ENDIF

            ! Type of friction law for each node
            ! ----------------------------------
            DO I=1, CF%DIM1
              NKFROT%I(I) = KFROT
            ENDDO
!
         ENDIF
!
!     Computation
!     -----------
!
      ELSE
!
         ! Friction by zones
         ! -----------------
         IF (FRICTB) THEN
! FH : FOR QUASI-BUBBLE
! FH : 2004/03/01
!jaj for quadratic elements
! =>
            IF (CF%ELM /= H%ELM) THEN
              IF (CF%ELM==12 .AND. H%ELM==11) THEN
                CALL FRICTION_BUBBLE
     &              (IKLE, NPOIN, NELEM, NELMAX, LINDNER, NKFROT,
     &               CHESTR, NDEFMA, LINDDP, LINDSP)
              ELSE IF (CF%ELM==13 .AND. H%ELM==11) THEN
                CALL FRICTION_QUAD
     &              (IKLE%I, NPOIN, NELEM, NELMAX, LINDNER, NKFROT,
     &               CHESTR, NDEFMA, LINDDP, LINDSP)
!                WRITE (LU,*) 
!     &           'FRICTION_CHOICE::QUADRATIC ELEMENTS NOT IMPLEMENTED.'
!                CALL PLANTE(1)
              ELSE
                WRITE (LU,*)
     &           'FRICTION_CHOICE::DISCRETISATION NOT IMPLEMENTED.'
                WRITE (LU,*)
     &           'CF%ELM, H%ELM: ',CF%ELM, H%ELM
                CALL PLANTE(1)
              ENDIF
            ENDIF
!
            CALL FRICTION_ZONES
     &           (MESH, H, U, V, S, CHESTR, CHBORD, NKFROT, NDEFMA,
     &            LINDDP, LINDSP, KFRO_B, NDEF_B, ITURB, LISRUG,
     &            LINDNER, VK, KARMAN, GRAV, T1, T2, CF, CFBOR)
! <=
!jaj for quadratic elements
! FH : 2004/03/01
! FH : FOR QUASI-BUBBLE
!
         ! Uniform friction
         ! ----------------
         ELSE
            CALL FRICTION_UNIF
     &           (MESH, H, U, V, CHESTR, S, KFROT, ITURB, LISRUG,
     &            LINDNER, SB, NDEF, DP, SP, VK, KARMAN, GRAV, T1,
     &            T2, CHBORD, CF, CFBOR)
!
         ENDIF   
      ENDIF
!
!======================================================================!
!======================================================================!
!
      RETURN
      END
