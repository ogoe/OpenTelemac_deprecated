!                    ************************
                     SUBROUTINE WALL_FRICTION
!                    ************************
!
     &(UETUTA,AUBOR,CFBOR,DISBOR,UN,VN,LIMPRO,NBOR,NPTFR,
     & KARMAN,PROPNU,LISRUG,KNEU,KDIR,KENT,KENTU,KADH,KLOG,IELMU,KP1BOR)
!
!***********************************************************************
! TELEMAC2D   V6P1                                   21/08/2010
!***********************************************************************
!
!brief    COMPUTES AUBOR, FRICTION ON BOUNDARIES.
!
!history  J-M HERVOUET (LNH)
!+        27/11/1992
!+
!+
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
!| AUBOR          |<--| LAW OF FRICTION ON BOUNDARIES
!|                |   | NUT*DU/DN=AUBOR*U+BUBOR
!| CFBOR          |<--| ADIMENSIONAL FRICTION COEFFICIENT ON BOUNDARIES
!| DISBOR         |-->| DISTANCE BETWEEN BOUNDARY AND NEAREST POINTS INSIDE
!| KADH           |-->| CONVENTION FOR NO SLIP BOUNDARY CONDITION
!| KARMAN         |-->| VON KARMAN CONSTANT
!| KDIR           |-->| CONVENTION FOR DIRICHLET POINT
!| KENT           |-->| CONVENTION FOR LIQUID INPUT WITH PRESCRIBED VALUE
!| KENTU          |-->| CONVENTION FOR LIQUID INPUT WITH PRESCRIBED VELOCITY
!| KLOG           |-->| CONVENTION FOR SOLID BOUNDARY
!| KNEU           |-->| CONVENTION FOR NEUMANN CONDITION
!| LIMPRO         |-->| BOUNDARY CONDITIONS FOR PROPAGATION (SEE PROPIN)
!| LISRUG         |-->| TURBULENCE REGIME (1: SMOOTH 2: ROUGH)
!| NBOR           |-->| GLOBAL NUMBER OF BOUNDARY POINTS
!| NPTFR          |-->| NUMBER OF BOUNDARY POINTS
!| PROPNU         |-->| LAMINAR DIFFUSION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN) :: NPTFR,LISRUG,KNEU,KDIR,KENT,KADH,KLOG,KENTU
      INTEGER, INTENT(IN) :: IELMU
      INTEGER, INTENT(IN) :: LIMPRO(NPTFR,6),NBOR(NPTFR),KP1BOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: CFBOR(*),UN(*),VN(*),DISBOR(*)
      DOUBLE PRECISION, INTENT(INOUT) :: AUBOR(*),UETUTA(*)
      DOUBLE PRECISION, INTENT(IN)    :: KARMAN,PROPNU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER N,K,IT
!
      DOUBLE PRECISION UTANG,DIST
!
!-----------------------------------------------------------------------
!
      INTRINSIC SQRT,LOG
!
!=======================================================================
!
!  LOOP ON THE BOUNDARY NODES
!
!  COMPUTES THE FRICTION VELOCITY ON THE WALL, AND AUBOR
!
!=======================================================================
!
!     ********************
      IF(LISRUG.EQ.1) THEN
!     ********************
!
        DO K=1,NPTFR
          N     = NBOR(K)
          UTANG = SQRT(UN(N)**2+VN(N)**2)
!         BEWARE : MODIFIED FROM PRINCIPLE NOTE
!         DIST  = DISBOR(K)*0.1D0
          DIST  = DISBOR(K)*0.33D0
!         COMPUTES UETOIL FOR SOLID BOUNDARIES
!         ----------------------------------------
!
!                            UETOIL
!         UETUTA REPRESENTS  ------ EVERYWHERE
!                            UTANG
!
!         UETUTA HAS THE ADVANTAGE OF BEING MEANINGFUL EVEN IF UTANG=0
!
!         INITIAL GUESS; THEN 5 ITERATIONS
          UETUTA(K) = 6.D-2
          DO IT=1,5
            IF(DIST*UETUTA(K)*UTANG/PROPNU .LT. 30.D0) THEN
              UETUTA(K) = 7.25D-2
            ELSE
              UETUTA(K)=1.D0/
     &                (5.5D0+LOG(DIST*UETUTA(K)*UTANG/PROPNU)/KARMAN)
            ENDIF
          ENDDO
        ENDDO
!
!     ************************
      ELSEIF(LISRUG.EQ.2) THEN
!     ************************
!
        DO K=1,NPTFR
          UETUTA(K) = SQRT( 0.5D0 * CFBOR(K) )
        ENDDO
!
!     ****
      ELSE
!     ****
!
        IF(LNG.EQ.1) WRITE(LU,400) LISRUG
        IF(LNG.EQ.2) WRITE(LU,401) LISRUG
400     FORMAT(1X,'WALL_FRICTION : REGIME DE TURBULENCE INCONNU :',1I6)
401     FORMAT(1X,'WALL_FRICTION: UNKNOWN TURBULENCE MODEL : ',1I6)
        CALL PLANTE(1)
        STOP
!
!     *****
      ENDIF
!     *****
!
! COMPUTES AUBOR
! --------------
!
!  AUBOR COUNTS FOR THE SEGMENT BETWEEN K AND KP1BOR(K)
!
!  LAW        : NUT * DU/DN = UETOIL**2 = -AUBOR*U(N+1)
!  CHANGED TO : NUT * DU/DN = UETOIL**2  *  U(N+1) / U(N)
!                           = UETOIL * (UETOIL/UTANG) * U(N+1)
!
      DO K=1,NPTFR
        IF(LIMPRO(K,5).EQ.KNEU) THEN
          N     = NBOR(K)
          UTANG = SQRT(UN(N)**2+VN(N)**2)
          AUBOR(K) = - UTANG * UETUTA(K)**2
        ELSE
          AUBOR(K) = 0.D0
        ENDIF
      ENDDO
!
!-----------------------------------------------------------------------
!
!     QUADRATIC ELEMENT HAS EXTRA POINTS ON THE BOUNDARY
!
      IF(IELMU.EQ.13) THEN
        DO K=1,NPTFR
          AUBOR(K+NPTFR) = (AUBOR(K)+AUBOR(KP1BOR(K)))*0.5D0
        ENDDO
      ENDIF
!
!=======================================================================
!
!                   /* END OF LOOP ON BOUNDARY NODES */
!
!=======================================================================
!
      RETURN
      END
