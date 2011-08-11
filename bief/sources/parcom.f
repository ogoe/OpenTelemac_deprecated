!                    *****************
                     SUBROUTINE PARCOM
!                    *****************
!
     &( X , ICOM , MESH )
!
!***********************************************************************
! BIEF   V6P1                                   21/08/2010
!***********************************************************************
!
!brief    COMPLEMENTS A VECTOR AT THE INTERFACES BETWEEN
!+                SUB-DOMAINS.
!+
!+            X CAN BE A BLOCK OF VECTORS. IN THIS CASE, ALL THE
!+                VECTORS IN THE BLOCK ARE TREATED.
!
!note     IMPORTANT : FROM RELEASE 5.9 ON, IDENTICAL VALUES ARE
!+                     ENSURED AT INTERFACE POINTS SO THAT DIFFERENT
!+                     PROCESSORS WILL ALWAYS MAKE THE SAME DECISION
!+                     IN TESTS ON REAL NUMBERS.
!
!warning  IF THE VECTORS HAVE A SECOND DIMENSION, IT IS
!+            IGNORED FOR THE TIME BEING
!
!history  J-M HERVOUET (LNHE)
!+        24/10/08
!+        V5P9
!+   AFTER REINHARD HINKELMANN (HANNOVER UNI.)
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        13/07/2010
!+        V6P0
!+   Translation of French comments within the FORTRAN sources into
!+   English comments
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        21/08/2010
!+        V6P0
!+   Creation of DOXYGEN tags for automated documentation and
!+   cross-referencing of the FORTRAN sources
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| ICOM           |-->| COMMUNICATION MODE
!|                |   | = 1 : VALUE WITH MAXIMUM ABSOLUTE VALUE TAKEN
!|                |   | = 2 : CONTRIBUTIONS ADDED
!|                |   | = 3 : MAXIMUM CONTRIBUTION RETAINED
!|                |   | = 4 : MINIMUM CONTRIBUTION RETAINED
!| MESH           |-->| MESH STRUCTURE
!| X              |<->| VECTOR OR BLOCK OF VECTORS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF, EX_PARCOM => PARCOM
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!-----------------------------------------------------------------------
!
      INTEGER, INTENT(IN) :: ICOM
!
!     STRUCTURES: VECTORS OR BLOCKS
!
      TYPE(BIEF_MESH), INTENT(INOUT)   :: MESH
      TYPE(BIEF_OBJ), INTENT(INOUT) :: X
!
!-----------------------------------------------------------------------
!
      TYPE(BIEF_OBJ), POINTER  :: X2,X3
      INTEGER NPOIN,NPLAN,IAN,NP11,NSEG
!
!***********************************************************************
!
!  OF NO USE IF A SUB-DOMAIN IS DISCONNECTED FROM THE OTHERS
!
      IF(NPTIR.EQ.0) RETURN
!
!-----------------------------------------------------------------------
!
      NPOIN = MESH%NPOIN
      NPLAN = 1
      IF(MESH%DIM.EQ.3) THEN
        NPOIN = BIEF_NBPTS(11,MESH)
        NPLAN = MESH%NPOIN/NPOIN
      ENDIF
!
!-----------------------------------------------------------------------
!
      IF(X%TYPE.EQ.2) THEN
!
!     VECTOR STRUCTURE
!
      IAN = 1
      CALL PARCOM2(X%R,X%R,X%R,NPOIN,NPLAN,ICOM,IAN,MESH)
!
      IF(X%ELM.EQ.13) THEN
        NP11=BIEF_NBPTS(11,MESH)
        NSEG=MESH%NSEG
        CALL PARCOM2_SEG(X%R(NP11+1:NP11+NSEG),
     &                   X%R(NP11+1:NP11+NSEG),
     &                   X%R(NP11+1:NP11+NSEG),
     &                   NSEG,1,ICOM,IAN,MESH,1)
      ENDIF
!
      ELSEIF(X%TYPE.EQ.4) THEN
!
!     BLOCK STRUCTURE
!
!     BEWARE: NUMBER LIMITED TO 3 |||||||||||||||||||||||||
      IAN = X%N
      IF(IAN.EQ.1) THEN
        X2 => X%ADR(1)%P
        X3 => X%ADR(1)%P
      ELSEIF(IAN.EQ.2) THEN
        X2 => X%ADR(2)%P
        X3 => X%ADR(2)%P
      ELSEIF(IAN.EQ.3) THEN
        X2 => X%ADR(2)%P
        X3 => X%ADR(3)%P
      ELSE
        IF(LNG.EQ.1) WRITE(LU,*) 'PARCOM PREVU JUSQU''A 3 VECTEURS'
        IF(LNG.EQ.2) WRITE(LU,*) 'PARCOM: NO MORE THAN 3 VECTORS'
        CALL PLANTE(1)
        STOP
      ENDIF
!
      CALL PARCOM2(X%ADR(1)%P%R,X2%R,X3%R,NPOIN,NPLAN,ICOM,IAN,MESH)
!
!     PROVISIONNALY 1 BY 1, COULD BE OPTIMISED
!
      IF(X%ADR(1)%P%ELM.EQ.13) THEN
        NP11=BIEF_NBPTS(11,MESH)
        NSEG=MESH%NSEG
        CALL PARCOM2_SEG(X%ADR(1)%P%R(NP11+1:NP11+NSEG),
     &                   X%ADR(1)%P%R(NP11+1:NP11+NSEG),
     &                   X%ADR(1)%P%R(NP11+1:NP11+NSEG),
!    *                   NSEG,1,ICOM,IAN,MESH)
     &                   NSEG,1,ICOM,1  ,MESH,1)
      ENDIF
      IF(IAN.GE.2.AND.X2%ELM.EQ.13) THEN
        NP11=BIEF_NBPTS(11,MESH)
        NSEG=MESH%NSEG
        CALL PARCOM2_SEG(X2%R(NP11+1:NP11+NSEG),
     &                   X2%R(NP11+1:NP11+NSEG),
     &                   X2%R(NP11+1:NP11+NSEG),
!    *                   NSEG,1,ICOM,IAN,MESH)
     &                   NSEG,1,ICOM,1  ,MESH,1)
      ENDIF
      IF(IAN.EQ.3.AND.X3%ELM.EQ.13) THEN
        NP11=BIEF_NBPTS(11,MESH)
        NSEG=MESH%NSEG
        CALL PARCOM2_SEG(X3%R(NP11+1:NP11+NSEG),
     &                   X3%R(NP11+1:NP11+NSEG),
     &                   X3%R(NP11+1:NP11+NSEG),
!    *                   NSEG,1,ICOM,IAN,MESH)
     &                   NSEG,1,ICOM,1  ,MESH,1)
      ENDIF
!
      ELSE
!
!     ERROR ON THE STRUCTURE
!
      IF (LNG.EQ.1) WRITE(LU,50) X%NAME,X%TYPE
      IF (LNG.EQ.1) WRITE(LU,53)
50    FORMAT(1X,'PARCOM (BIEF) : NOM DE X : ',A6,'  TYPE : ',1I6)
53    FORMAT(1X,'                CAS NON PREVU')
      IF (LNG.EQ.2) WRITE(LU,51) X%NAME,X%TYPE
      IF (LNG.EQ.2) WRITE(LU,54)
51    FORMAT(1X,'PARCOM (BIEF): NAME OF X: ',A6,'  TYPE : ',1I6)
54    FORMAT(1X,'               UNEXPECTED CASE')
      CALL PLANTE(1)
      STOP
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      RETURN
      END
