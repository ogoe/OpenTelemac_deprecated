!                  *****************************
                   SUBROUTINE SUSPENSION_EROSION
!                  *****************************
!
     &(TAUP,HN,FDM,AVA,NPOIN,CHARR,XMVE,XMVS,VCE,GRAV,HMIN,XWC,
     & ZERO,ZREF,AC,FLUER,CSTAEQ,QSC,ICQ,DEBUG)
!
!***********************************************************************
! SISYPHE   V6P1                                   21/07/2011
!***********************************************************************
!
!brief    COMPUTES THE FLUX OF DEPOSITION AND EROSION.
!
!history  J-M HERVOUET + C VILLARET
!+        17/09/2009
!+
!+
!
!history  CV
!+        04/05/2010
!+        V6P0
!+   MODIFICATION FOR FREDSOE: EQUILIBRIUM CONCENTRATIONS
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
!history  C.VILLARET (EDF-LNHE), P.TASSI (EDF-LNHE)
!+        19/07/2011
!+        V6P1
!+   Name of variables   
!+   
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| AC             |<->| CRITICAL SHIELDS PARAMETER (SHOULD BE INPUT ONLY)
!| ACLADM         |-->| MEAN DIAMETER OF SEDIMENT
!| AVA            |-->| VOLUME PERCENT OF SEDIMENT CLASS IN THE TOP ACTIVE LAYER
!| CHARR          |-->| LOGICAL, IF BEDLOAD OR NOT
!| CSTAEQ         |<->| EQUILIBRIUM CONCENTRATION
!| DEBUG          |-->| FLAG FOR DEBUGGING
!| FLUER          |<->| EROSION FLUX
!| GRAV           |-->| ACCELERATION OF GRAVITY
!| HMIN           |-->| MINIMUM VALUE OF WATER DEPTH
!| HN             |-->| WATER DEPTH
!| ICQ            |-->| REFERENCE CONCENTRATION FORMULA
!| NPOIN          |-->| NUMBER OF POINTS
!| QSC            |-->| BED LOAD TRANSPORT RATE
!| TAUP           |-->| SKIN FRICTION
!| VCE            |-->| FLOW VISCOSITY
!| XMVE           |-->| FLUID DENSITY 
!| XMVS           |-->| WATER DENSITY
!| XWC            |-->| SETTLING VELOCITIES 
!| ZERO           |-->| ZERO
!| ZREF           |-->| REFERENCE ELEVATION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE INTERFACE_SISYPHE, EX_SUSPENSION_EROSION=>SUSPENSION_EROSION
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU

      ! 2/ GLOBAL VARIABLES
      ! -------------------
      TYPE (BIEF_OBJ),  INTENT(IN)    :: TAUP,HN,ZREF,QSC
      INTEGER,          INTENT(IN)    :: NPOIN,DEBUG,ICQ
      LOGICAL,          INTENT(IN)    :: CHARR
      DOUBLE PRECISION, INTENT(IN)    :: XMVE,XMVS,GRAV,HMIN,XWC,VCE
      DOUBLE PRECISION, INTENT(IN)    :: AVA(NPOIN),ZERO,FDM
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: FLUER,CSTAEQ
      DOUBLE PRECISION, INTENT(INOUT) :: AC

      ! 3/ LOCAL VARIABLES
      ! -------------------

      INTEGER I
!
!======================================================================!
!======================================================================!
!                               PROGRAM                                !
!======================================================================!
!======================================================================!
!
!  COMPUTES THE NEAR BED EQUILIBRIUM CONCENTRATION --> CSTAEQ (MEAN DIAMETER)
!
      IF(ICQ.EQ.1) THEN
!
        IF(DEBUG > 0) WRITE(LU,*) 'SUSPENSION_FREDSOE'
        CALL SUSPENSION_FREDSOE(FDM,TAUP,NPOIN,
     &                          GRAV,XMVE,XMVS,ZERO,AC,CSTAEQ)
        IF(DEBUG > 0) WRITE(LU,*) 'END SUSPENSION_FREDSOE'
!
        DO I=1,NPOIN
          CSTAEQ%R(I)=CSTAEQ%R(I)*AVA(I)
        ENDDO
        CALL OS('X=CY    ', X=FLUER, Y=CSTAEQ, C=XWC)
!
      ELSEIF(ICQ.EQ.2) THEN
!
        IF(DEBUG > 0) WRITE(LU,*) 'SUSPENSION_BIJKER'
        CALL SUSPENSION_BIJKER(TAUP,HN,NPOIN,CHARR,QSC,ZREF,
     &                         ZERO,HMIN,CSTAEQ,XMVE)
        IF(DEBUG > 0) WRITE(LU,*) 'END SUSPENSION_BIJKER'
!       NO MULTIPLICATION BY AVA BECAUSE AVA HAS ALREADY BEEN TAKEN
!       INTO ACCOUNT IN THE BEDLOAD TRANSPORT RATE
        CALL OS('X=CY    ', X=FLUER, Y=CSTAEQ, C=XWC)
!
      ELSEIF(ICQ.EQ.3) THEN
!
        IF(DEBUG > 0) WRITE(LU,*) 'SUSPENSION_VANRIJN' 
        CALL SUSPENSION_VANRIJN(FDM,TAUP,NPOIN,GRAV,XMVE,XMVS,
     &                          VCE,ZERO,AC,CSTAEQ,ZREF)
        IF(DEBUG > 0) WRITE(LU,*) 'END SUSPENSION_VANRIJN'
        DO I=1,NPOIN
          CSTAEQ%R(I)=CSTAEQ%R(I)*AVA(I)
        ENDDO
        CALL OS('X=CY    ', X=FLUER, Y=CSTAEQ, C=XWC) 
!
      ENDIF 
!
!======================================================================!
!======================================================================!
!
      RETURN
      END
C
C#######################################################################
C
