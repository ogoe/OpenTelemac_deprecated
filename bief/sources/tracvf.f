!                    *****************
                     SUBROUTINE TRACVF
!                    *****************
!
     &(F,FN,FSCEXP,H,HN,FXMAT,FXMATPAR,
     & V2DPAR,UNSV2D,DDT,FXBOR,FBOR,SMH,YASMH,T1,T2,T4,T5,T6,T7,T8,
     & MESH,LIMTRA,KDIR,KDDL,OPTSOU,IOPT2,FLBORTRA,MSK,DT)
!
!***********************************************************************
! BIEF   V6P1                                   21/08/2010
!***********************************************************************
!
!brief    COMPUTES THE TRACER FOR FINITE VOLUME SCHEME.
!+                TO COMPLETE.
!
!history  C-T PHAM (LNHE)
!+        06/02/09
!+        V6P0
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
!| DDT            |-->| SUB TIME-STEP
!| DT             |-->| TIME-STEP
!| F              |<--| VALUES OF F AT TIME N+1.
!| FBOR           |-->| VALUES OF F AT THE PRESCRIBED BOUNDARIES
!| FLBORTRA       |<->| FLUX OF TRACER AT THE BOUNDARIES
!| FN             |-->| VALEURS DU TRACEUR A L'ETAPE N.
!| FSCEXP         |-->| EXPLICIT SOURCE TERM FOR F
!| FXBOR          |-->| FLUXES ON BOUNDARIES
!| FXMAT          |-->| FLUXES (NON ASSEMBLED IN PARALLEL)
!| FXMATPAR       |-->| FLUXES (ASSEMBLED IN PARALLEL)
!| H              |-->| WATER DEPTH AT TIME N+1
!| HN             |-->| WATER DEPTH AT TIME N
!| IOPT2          |-->| 0 : UCONV OBEYS CONTINUITY EQUATION
!|                |   | 1 : UCONV DOES NOT OBEY CONTINUITY EQUATION
!| KDDL           |-->| CONVENTION FOR DEGREE OF FREEDOM
!| KDIR           |-->| CONVENTION FOR DIRICHLET POINT
!| LIMTRA         |-->| TECHNICAL BOUNDARY CONDITIONS FOR TRACERS
!| MESH           |-->| MESH STRUCTURE
!| MSK            |-->| IF YES, MASKING OF DRY ELEMENTS
!| OPTSOU         |-->| TYPE OF SOURCES
!|                |   | 1: NORMAL
!|                |   | 2: DIRAC
!| SMH            |-->| SOURCE TERM IN CONTINUITY EQUATION
!| T1             |<->| BIEF_OBJ STRUCTURE USED AS WORK ARRAY
!| T2             |<->| BIEF_OBJ STRUCTURE USED AS WORK ARRAY
!| T4             |<->| BIEF_OBJ STRUCTURE USED AS WORK ARRAY
!| T5             |<->| BIEF_OBJ STRUCTURE USED AS WORK ARRAY
!| T6             |<->| BIEF_OBJ STRUCTURE USED AS WORK ARRAY
!| T7             |<->| BIEF_OBJ STRUCTURE USED AS WORK ARRAY
!| T8             |<->| BIEF_OBJ STRUCTURE USED AS WORK ARRAY
!| UNSV2D         |-->| INVERSE OF INTEGRALS OF TEST FUNCTIONS
!| V2DPAR         |-->| INTEGRAL OF TEST FUNCTIONS, ASSEMBLED IN PARALLEL
!| YASMH          |-->| IF YES, SMH MUST BE TAKEN INTO ACCOUNT
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF, EX_TRACVF => TRACVF
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN)           :: KDIR,KDDL,OPTSOU,LIMTRA(*)
      INTEGER, INTENT(IN)           :: IOPT2
      DOUBLE PRECISION, INTENT(IN)  :: DDT,DT
      TYPE(BIEF_OBJ), INTENT(INOUT) :: F,T1,T2,T4,T5,T6,T7,T8,FLBORTRA
      TYPE(BIEF_OBJ), INTENT(IN)    :: FN,H,HN,V2DPAR,SMH,FBOR,FSCEXP
      TYPE(BIEF_OBJ), INTENT(IN)    :: FXBOR,UNSV2D
      DOUBLE PRECISION, INTENT(IN)  :: FXMAT(*),FXMATPAR(*)
      TYPE(BIEF_MESH), INTENT(INOUT):: MESH
      LOGICAL, INTENT(IN)           :: YASMH,MSK
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER I
!
!-----------------------------------------------------------------------
!
      IF(IOPT2.EQ.0) THEN
!
!-----------------------------------------------------------------------
!
!     CASE WHERE THE ADVECTION FIELD SATISFIES THE CONTINUITY EQUATION
!     (THE DEPTH COULD BE COMPUTED BY INTERPOLATION IN TIME)
!
!     T4 WILL TAKE THE SUCCESSIVE VALUES OF F (INITIALISED IN CVTRVF)
!
      CALL TVF(F%R,FN%R,T4%R,T5%R,FXMAT,FXMATPAR,UNSV2D%R,DDT,
     &         FXBOR%R,T7%R,T8,FBOR%R,SMH%R,YASMH,FSCEXP%R,
     &         MESH%NSEG,MESH%NPOIN,MESH%NPTFR,
     &         MESH%GLOSEG%I,MESH%GLOSEG%DIM1,
     &         MESH%NBOR%I,LIMTRA,KDIR,KDDL,
     &         OPTSOU,T5%R,IOPT2,FLBORTRA%R,DDT/DT,MESH,F)
!
!-----------------------------------------------------------------------
!
!     CASE WHERE THE ADVECTION FIELD DOES NOT SATISFY THE CONTINUITY EQUATION
!
      ELSEIF(IOPT2.EQ.1) THEN
!
!     T1 WILL TAKE THE SUCCESSIVE VALUES OF HN COMPUTED WITH CONTINUITY
!     T2 WILL TAKE THE SUCCESSIVE VALUES OF H COMPUTED WITH CONTINUITY
!     T4 WILL TAKE THE SUCCESSIVE VALUES OF F
!     T5 WILL TAKE THE SUCCESSIVE VALUES OF TRUE DEPTH
!
!     H2 DEPTH BY CONTINUITY EQUATION
!
      CALL HVF(T2%R,T1%R,FXMAT,UNSV2D%R,DDT,T7%R,SMH%R,
     &         YASMH,MESH%NSEG,MESH%NPOIN,MESH%NPTFR,
     &         MESH%GLOSEG%I,MESH%GLOSEG%DIM1,MESH%NBOR%I,OPTSOU,
     &         T8,MESH,MSK)
!
      CALL TVF(F%R,FN%R,T4%R,T2%R,FXMAT,FXMATPAR,UNSV2D%R,DDT,
     &         FXBOR%R,T7%R,T8,FBOR%R,SMH%R,YASMH,FSCEXP%R,
     &         MESH%NSEG,MESH%NPOIN,MESH%NPTFR,
     &         MESH%GLOSEG%I,MESH%GLOSEG%DIM1,
     &         MESH%NBOR%I,LIMTRA,KDIR,KDDL,
     &         OPTSOU,T5%R,IOPT2,FLBORTRA%R,DDT/DT,MESH,F)
!
!-----------------------------------------------------------------------
!
      ELSE
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'TRACVF : OPTION INCONNUE'
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) 'TRACVF: UNKNOWN OPTION'
        ENDIF
        CALL PLANTE(1)
        STOP
      ENDIF
!
!-----------------------------------------------------------------------
!
      RETURN
      END
