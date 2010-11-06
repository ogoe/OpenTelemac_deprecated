C                       **************************
                        SUBROUTINE FRICTION_BUBBLE
C                       **************************
C
     & (IKLE, NPOIN, NELEM, NELMAX, LINDNER, NKFROT, CHESTR, NDEFMA,
     &  LINDDP, LINDSP)
C
C***********************************************************************
C  TELEMAC-2D VERSION 5.5                 J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C 22/12/04 : subroutine written by F. Huvelin
C
C
C   ! --------------------------------------------------------------- !
C   ! Computation of the friction vector for the quasi-bubble element !
C   ! --------------------------------------------------------------- !
C
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

!=======================================================================!
!=======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS                !
!=======================================================================!
!=======================================================================!

      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ), INTENT(IN)    :: IKLE
      INTEGER,        INTENT(IN)    :: NPOIN, NELEM, NELMAX
      LOGICAL,        INTENT(IN)    :: LINDNER
      TYPE(BIEF_OBJ), INTENT(INOUT) :: NKFROT, CHESTR, NDEFMA
      TYPE(BIEF_OBJ), INTENT(INOUT) :: LINDDP, LINDSP
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER :: I, I1, I2, I3
!
!=======================================================================!
!=======================================================================!
!                               PROGRAMME                               !
!=======================================================================!
!=======================================================================!
!
      DO I = NPOIN + 1, NPOIN + NELEM
!
         I1 = IKLE%I(I - NPOIN           )
         I2 = IKLE%I(I - NPOIN +   NELMAX)
         I3 = IKLE%I(I - NPOIN + 2*NELMAX)
!
         ! COMPUTING THE VALUES OF THE MIDDLE-NODE
         ! ---------------------------------------
         IF (NKFROT%I(I1).EQ.NKFROT%I(I2)) THEN

            ! The 3 nodes have the same law !
            ! ***************************** !
            IF (NKFROT%I(I1).EQ.NKFROT%I(I3))THEN
!
               NKFROT%I(I) = NKFROT%I(I1)
               CHESTR%R(I) = (CHESTR%R(I3) + CHESTR%R(I2) +CHESTR%R(I1))
     &                     / 3.D0
               NDEFMA%R(I) = (NDEFMA%R(I3) + NDEFMA%R(I2) +NDEFMA%R(I1))
     &                     / 3.D0
!
               IF (LINDNER) THEN
                  LINDDP%R(I) = ( LINDDP%R(I3) + LINDDP%R(I2)
     &                           +LINDDP%R(I1) )/3.D0
!
                  LINDSP%R(I) = ( LINDSP%R(I3) + LINDSP%R(I2)
     &                           +LINDSP%R(I1) )/3.D0
               ENDIF
!
            ! The nodes "1" and "2" have the same law !
            ! *************************************** !
            ELSE
               NKFROT%I(I) = NKFROT%I(I1)
               CHESTR%R(I) = (CHESTR%R(I2) + CHESTR%R(I1))/2.D0
               NDEFMA%R(I) = (NDEFMA%R(I2) + NDEFMA%R(I1))/2.D0
!
               IF (LINDNER) THEN
                  LINDDP%R(I) = (LINDDP%R(I2) + LINDDP%R(I1))/2.D0
                  LINDSP%R(I) = (LINDSP%R(I2) + LINDSP%R(I1))/2.D0
               ENDIF
            ENDIF
!
         ! The nodes "2" and "3" have the same law !
         ! *************************************** !
         ELSE IF (NKFROT%I(I2).EQ.NKFROT%I(I3)) THEN
!
            NKFROT%I(I) = NKFROT%I(I2)
            CHESTR%R(I) = (CHESTR%R(I3) + CHESTR%R(I2))/2.D0
            NDEFMA%R(I) = (NDEFMA%R(I3) + NDEFMA%R(I2))/2.D0
!
            IF (LINDNER) THEN
               LINDDP%R(I) = (LINDDP%R(I3) + LINDDP%R(I2))/2.D0
               LINDSP%R(I) = (LINDSP%R(I3) + LINDSP%R(I2))/2.D0
            ENDIF
!
         ! The 3 nodes have different laws : value of the node "1" kept !
         ! ************************************************************ !
         ELSE
            NKFROT%I(I) = NKFROT%I(I1)
            CHESTR%R(I) = CHESTR%R(I1)
            NDEFMA%R(I) = NDEFMA%R(I1)
!
            IF (LINDNER) THEN
               LINDDP%R(I) = LINDDP%R(I1)
               LINDSP%R(I) = LINDSP%R(I1)
            ENDIF
         ENDIF
      ENDDO
!
!=======================================================================!
!=======================================================================!
!
      RETURN
      END
