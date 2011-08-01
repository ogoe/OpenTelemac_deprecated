C                       ************************
                        SUBROUTINE FRICTION_QUAD
C                       ************************
C
     & (IKLE, NPOIN, NELEM, NELMAX, LINDNER, NKFROT, CHESTR, NDEFMA,
     &  LINDDP, LINDSP)
C
C***********************************************************************
C  TELEMAC-2D VERSION 6.0                          JACEK JANKOWSKI (BAW)
C***********************************************************************
C
C 22/12/04 : subroutine written by F. Huvelin
C Fri Apr 17 12:58:01 CEST 2009 jaj pinxit 
C
C
C   ! --------------------------------------------------------------- !
C   ! Computation of the friction vector for the quadratic element    !
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
      INTEGER,        INTENT(IN)    :: NPOIN,NELEM,NELMAX
      INTEGER,        INTENT(IN)    :: IKLE(NELMAX,6)
      LOGICAL,        INTENT(IN)    :: LINDNER
      TYPE(BIEF_OBJ), INTENT(INOUT) :: NKFROT,CHESTR,NDEFMA
      TYPE(BIEF_OBJ), INTENT(INOUT) :: LINDDP,LINDSP
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER :: IELEM
!
!=======================================================================!
!=======================================================================!
!                               PROGRAMME                               !
!=======================================================================!
!=======================================================================!
!
! if the 11 nodes have the same friction law, interpolation a la cg1113
! for the values at 13 additional nodes in the middle of the edges 
!
!        X(IKLE(IELEM,4)) = 0.5D0 * ( X(IKLE(IELEM,1))
!     &                             + X(IKLE(IELEM,2)) )
!        X(IKLE(IELEM,5)) = 0.5D0 * ( X(IKLE(IELEM,2))
!     &                             + X(IKLE(IELEM,3)) )
!        X(IKLE(IELEM,6)) = 0.5D0 * ( X(IKLE(IELEM,3)) 
!     &                             + X(IKLE(IELEM,1)) )
!
! well, if the the friction laws differ, take the value on the previous 
! node by circumventing the element... 
!
!        X(IKLE(IELEM,4)) = X(IKLE(IELEM,1))
!        X(IKLE(IELEM,5)) = X(IKLE(IELEM,2))
!        X(IKLE(IELEM,6)) = X(IKLE(IELEM,3)) 
!
! assumed the trivial case of -one- zone only is naturally excluded
!   
      DO IELEM = 1,NELEM

        IF(NKFROT%I(IKLE(IELEM,1)).EQ.NKFROT%I(IKLE(IELEM,2))) THEN
          NKFROT%I(IKLE(IELEM,4)) = NKFROT%I(IKLE(IELEM,1))
          CHESTR%R(IKLE(IELEM,4)) = 
     &      0.5D0*(CHESTR%R(IKLE(IELEM,1))+CHESTR%R(IKLE(IELEM,2)))
          NDEFMA%R(IKLE(IELEM,4)) = 
     &      0.5D0*(NDEFMA%R(IKLE(IELEM,1))+NDEFMA%R(IKLE(IELEM,2)))
        ELSE 
          NKFROT%I(IKLE(IELEM,4)) = NKFROT%I(IKLE(IELEM,1))
          CHESTR%R(IKLE(IELEM,4)) = CHESTR%R(IKLE(IELEM,1))
          NDEFMA%R(IKLE(IELEM,4)) = NDEFMA%R(IKLE(IELEM,1))
        ENDIF

        IF(NKFROT%I(IKLE(IELEM,2)).EQ.NKFROT%I(IKLE(IELEM,3))) THEN
          NKFROT%I(IKLE(IELEM,5)) = NKFROT%I(IKLE(IELEM,2))
          CHESTR%R(IKLE(IELEM,5)) = 
     &      0.5D0*(CHESTR%R(IKLE(IELEM,2))+CHESTR%R(IKLE(IELEM,3)))
          NDEFMA%R(IKLE(IELEM,5)) = 
     &      0.5D0*(NDEFMA%R(IKLE(IELEM,2))+NDEFMA%R(IKLE(IELEM,3)))
        ELSE 
          NKFROT%I(IKLE(IELEM,5)) = NKFROT%I(IKLE(IELEM,2))
          CHESTR%R(IKLE(IELEM,5)) = CHESTR%R(IKLE(IELEM,2))
          NDEFMA%R(IKLE(IELEM,5)) = NDEFMA%R(IKLE(IELEM,2))
        ENDIF

        IF(NKFROT%I(IKLE(IELEM,3)).EQ.NKFROT%I(IKLE(IELEM,1))) THEN
          NKFROT%I(IKLE(IELEM,6)) = NKFROT%I(IKLE(IELEM,3))
          CHESTR%R(IKLE(IELEM,6)) = 
     &      0.5D0*(CHESTR%R(IKLE(IELEM,3))+CHESTR%R(IKLE(IELEM,1)))
          NDEFMA%R(IKLE(IELEM,6)) = 
     &      0.5D0*(NDEFMA%R(IKLE(IELEM,3))+NDEFMA%R(IKLE(IELEM,1)))
        ELSE 
          NKFROT%I(IKLE(IELEM,6)) = NKFROT%I(IKLE(IELEM,3))
          CHESTR%R(IKLE(IELEM,6)) = CHESTR%R(IKLE(IELEM,3))
          NDEFMA%R(IKLE(IELEM,6)) = NDEFMA%R(IKLE(IELEM,3))
        ENDIF

      ENDDO

      ! this rare case separately 

      IF (LINDNER) THEN 

        DO IELEM = 1,NELEM

          IF (NKFROT%I(IKLE(IELEM,1)).EQ.NKFROT%I(IKLE(IELEM,2))) THEN
            LINDDP%R(IKLE(IELEM,4)) = 
     &        0.5D0*(LINDDP%R(IKLE(IELEM,1))+LINDDP%R(IKLE(IELEM,2)))
            LINDSP%R(IKLE(IELEM,4)) = 
     &        0.5D0*(LINDSP%I(IKLE(IELEM,1))+LINDSP%R(IKLE(IELEM,2)))
          ELSE 
            LINDDP%R(IKLE(IELEM,4)) = LINDDP%R(IKLE(IELEM,1))
            LINDSP%R(IKLE(IELEM,4)) = LINDSP%R(IKLE(IELEM,1))
          ENDIF

          IF (NKFROT%I(IKLE(IELEM,2)).EQ.NKFROT%I(IKLE(IELEM,3))) THEN
            LINDDP%R(IKLE(IELEM,5)) = 
     &        0.5D0*(LINDDP%R(IKLE(IELEM,2))+LINDDP%R(IKLE(IELEM,3)))
            LINDSP%R(IKLE(IELEM,5)) = 
     &        0.5D0*(LINDSP%I(IKLE(IELEM,2))+LINDSP%R(IKLE(IELEM,3)))
          ELSE 
            LINDDP%R(IKLE(IELEM,5)) = LINDDP%R(IKLE(IELEM,2))
            LINDSP%R(IKLE(IELEM,5)) = LINDSP%R(IKLE(IELEM,2))
          ENDIF

          IF (NKFROT%I(IKLE(IELEM,3)).EQ.NKFROT%I(IKLE(IELEM,1))) THEN
            LINDDP%R(IKLE(IELEM,6)) = 
     &        0.5D0*(LINDDP%R(IKLE(IELEM,3))+LINDDP%R(IKLE(IELEM,1)))
            LINDSP%R(IKLE(IELEM,6)) = 
     &        0.5D0*(LINDSP%I(IKLE(IELEM,3))+LINDSP%R(IKLE(IELEM,1)))
          ELSE 
            LINDDP%R(IKLE(IELEM,6)) = LINDDP%R(IKLE(IELEM,3))
            LINDSP%R(IKLE(IELEM,6)) = LINDSP%R(IKLE(IELEM,3))
          ENDIF

        ENDDO 

      ENDIF
!
!=======================================================================!
!=======================================================================!
!
      RETURN
      END SUBROUTINE FRICTION_QUAD
