!                    *****************
                     SUBROUTINE KEPSCL
!                    *****************
!
     &(KBOR,EBOR,AUBOR,CF,CFBOR,DISBOR,
     & UN,VN,HN,LIMKEP,LIUBOR,LIMPRO,NBOR,NPTFR,
     & KARMAN,CMU,C2,ESTAR,SCHMIT,LISRUG,PROPNU,KMIN,EMIN,
     & KNEU,KDIR,KENT,KENTU,KADH,KLOG,UETUTA)
!
!***********************************************************************
! TELEMAC2D   V6P1                                   21/08/2010
!***********************************************************************
!
!brief    COMPUTES KBOR, EBOR AND AUBOR WHEN THE TURBULENCE
!+                MODEL IS K-EPSILON.
!
!history  J-M HERVOUET (LNH)
!+        27/11/1992
!+
!+
!
!history  L. VAN HAREN (LNH)
!+        26/04/1994
!+        V5P2
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
!| CF             |-->| ADIMENSIONAL FRICTION COEFFICIENT
!| CFBOR          |-->| ADIMENSIONAL FRICTION COEFFICIENT OF BOUNDARIES
!| CMU            |-->| CONSTANT OF K-EPSILON MODEL
!| C2             |-->| CONSTANT OF K-EPSILON MODEL
!| DISBOR         |-->| DISTANCE TO BOUNDARY
!| EBOR           |<--| TURBULENT ENERGY DISSIPATION AT BOUNDARY
!| EMIN           |-->| MINIMUM TURBULENT ENERGY DISSIPATION
!| ESTAR          |-->| CONSTANT OF K-EPSILON MODEL
!| HN             |-->| WATER DEPTH AT TIME T(N)
!| KADH           |-->| CONVENTION FOR NO SLIP BOUNDARY CONDITION
!| KARMAN         |-->| VON KARMAN CONSTANT
!| KBOR           |<--| TURBULENTE KINETIC ENERGY ON BOUNDARIES
!| KDIR           |-->| CONVENTION FOR DIRICHLET POINT
!| KENT           |-->| CONVENTION FOR LIQUID INPUT WITH PRESCRIBED VALUE
!| KENTU          |-->| CONVENTION FOR LIQUID INPUT WITH PRESCRIBED VELOCITY
!| KLOG           |-->| CONVENTION FOR SOLID BOUNDARY
!| KMIN           |-->| MINIMUM TURBULENT KINETIC ENERGY
!| KNEU           |-->| CONVENTION FOR NEUMANN CONDITION
!| LIMKEP         |-->| BOUNDARY CONDITIONS ON K AND EPSILON
!| LIMPRO         |-->| BOUNDARY CONDITIONS FOR PROPAGATION
!| LISRUG         |-->| TURBULENCE REGIME (1: SMOOTH 2: ROUGH)
!| LIUBOR         |-->| TYPE OF BOUNDARY CONDITIONS ON VELOCITY
!| NBOR           |-->| GLOBAL NUMBER OF BOUNDARY POINTS
!| NPTFR          |-->| NUMBER OF BOUNDARY POINTS
!| PROPNU         |-->| MOLECULAR DIFFUSION
!| SCHMIT         |-->| CONSTANT OF K-EPSILON MODEL
!| UETUTA         |-->| FRICTION VELOCITY/TANGENTIAL VELOCITY
!| UN             |<->| X-COMPONENT OF VELOCITY AT TIME T(N)
!| VN             |<->| Y-COMPONENT OF VELOCITY AT TIME T(N)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN) :: NPTFR,LISRUG
      INTEGER, INTENT(IN) :: KNEU,KDIR,KENT,KADH,KLOG,KENTU
      INTEGER, INTENT(IN) :: LIMPRO(NPTFR,6),NBOR(NPTFR)
      INTEGER, INTENT(IN) :: LIMKEP(NPTFR,2),LIUBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: KMIN,EMIN
      DOUBLE PRECISION, INTENT(IN)    :: CF(*),CFBOR(*),UETUTA(*)
      DOUBLE PRECISION, INTENT(IN)    :: UN(*),VN(*),HN(*),DISBOR(*)
      DOUBLE PRECISION, INTENT(INOUT) :: KBOR(*),EBOR(*),AUBOR(*)
      DOUBLE PRECISION, INTENT(IN) :: KARMAN,CMU,C2,ESTAR,SCHMIT,PROPNU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER N,K
!
      DOUBLE PRECISION KFOND,EFOND,TIERS,CEPS,USTAR
      DOUBLE PRECISION SSQCMU,UTANG,DIST,DENOM,EBORD,KBORD
!
!-----------------------------------------------------------------------
!
      INTRINSIC SQRT,MAX
!
!-----------------------------------------------------------------------
!
      TIERS  = 1.D0/3.D0
      SSQCMU = 1.D0/ SQRT(CMU)
!
!=======================================================================
!
!  LOOP ON THE BOUNDARY NODES
!
!  COMPUTES KBOR,EBOR
!
!=======================================================================
!
      DO K=1,NPTFR
!
         KBOR(K) = 0.D0
         EBOR(K) = 0.D0
         N     = NBOR(K)
         UTANG = SQRT( UN(N)**2 + VN(N)**2 )
!        BEWARE : MODIFIED FROM PRINCIPLE NOTE
!        DIST  = DISBOR(K)*0.1D0
         DIST  = DISBOR(K)*0.33D0
!
!  DIRICHLET ON K
!  ---------------
!
         IF(LIMKEP(K,1).EQ.KDIR) THEN
!        ----------------------------
!
!           ************************************************
            IF(LIUBOR(K).EQ.KENT.OR.LIUBOR(K).EQ.KENTU) THEN
!           ************************************************
!
!              INPUT BOUNDARY: TURBULENCE DUE TO THE BOTTOM
!
               CEPS    = C2*SQRT(CMU)/SQRT(ESTAR*SCHMIT) /
     &                   (0.5D0*CF(N))**0.75D0
               DENOM   = CEPS * 0.5D0*CF(N)
               USTAR = SQRT( 0.5D0 * CF(N) * ( UN(N)**2 + VN(N)**2 ) )
               KBOR(K) = C2 * USTAR**2 / MAX(DENOM,1.D-10)
!
!           ***************************************************
            ELSEIF(LIUBOR(K).EQ.KLOG.OR.LIUBOR(K).EQ.KADH) THEN
!           ***************************************************
!
!              WALL
!
               CEPS    = C2*SQRT(CMU)/SQRT(ESTAR*SCHMIT) /
     &                   (0.5D0*CF(N))**0.75D0
               DENOM   = CEPS * 0.5D0*CF(N)
               USTAR = SQRT( 0.5D0 * CF(N) * ( UN(N)**2 + VN(N)**2 ) )
               KFOND   = C2 * USTAR**2 / MAX(DENOM,1.D-10)
               KBORD   = SSQCMU*(UETUTA(K)*UTANG)**2
               KBOR(K) = KBORD + KFOND
!
!           ****
            ELSE
!           ****
!
               IF(LNG.EQ.1) WRITE(LU,500) K,LIUBOR(K)
               IF(LNG.EQ.2) WRITE(LU,501) K,LIUBOR(K)
500            FORMAT(1X,'KEPSCL: POINT DE BORD ',1I6,
     &                   'CAS NON PREVU POUR KBOR',1X,'LIUBOR=',1I6)
501            FORMAT(1X,'KEPSCL: BOUNDARY POINT ',1I6,
     &                   'UNKNOWN CASE FOR KBOR',1X,'LIUBOR=',1I6)
               CALL PLANTE(1)
               STOP
!
!           *****
            ENDIF
!           *****
!
         ENDIF
!        -----
!
!  DIRICHLET ON EPSILON
!  ---------------------
!
         IF(LIMKEP(K,2).EQ.KDIR) THEN
!        ----------------------------
!
!           ************************************************
            IF(LIUBOR(K).EQ.KENT.OR.LIUBOR(K).EQ.KENTU) THEN
!           ************************************************
!
!              INPUT BOUNDARY: TURBULENCE DUE TO THE BOTTOM
!
               DENOM   = SQRT(0.5D0*CF(N)) * HN(N)
               USTAR   = SQRT(0.5D0*CF(N) * ( UN(N)**2 + VN(N)**2 ) )
               EFOND   = USTAR**3 / MAX(DENOM,1.D-10)
               EBOR(K) = MAX( EFOND , EMIN )
!
!           ***************************************************
            ELSEIF(LIUBOR(K).EQ.KLOG.OR.LIUBOR(K).EQ.KADH) THEN
!           ***************************************************
!
!              WALL
!
               DENOM   = SQRT(0.5D0*CF(N)) * HN(N)
               USTAR   = SQRT(0.5D0*CF(N) * ( UN(N)**2 + VN(N)**2 ) )
               EFOND   = USTAR**3 / MAX(DENOM,1.D-10)
               EBORD   = (UETUTA(K)*UTANG)**3 / ( KARMAN*DIST )
               EBOR(K) = MAX( EBORD + EFOND, EMIN )
!
!           ****
            ELSE
!           ****
!
!              OTHER
!
               IF(LNG.EQ.1) WRITE(LU,600) K,LIUBOR(K)
               IF(LNG.EQ.2) WRITE(LU,601) K,LIUBOR(K)
600            FORMAT(1X,'KEPSCL: POINT DE BORD ',1I6,
     &                   'CAS NON PREVU POUR EBOR',1X,'LIUBOR=',1I6)
601            FORMAT(1X,'KEPSCL: BOUNDARY POINT ',1I6,
     &                   'UNKNOWN CASE FOR EBOR',1X,'LIUBOR=',1I6)
               CALL PLANTE(1)
               STOP
!
!           *****
            ENDIF
!           *****
!
         ENDIF
!        -----
!
      ENDDO
!
!=======================================================================
!
!                   /* END OF LOOP ON BOUNDARY NODES */
!
!=======================================================================
!
      RETURN
      END
