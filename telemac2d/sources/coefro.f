! modif CV sept 2011
!                    *****************
                     SUBROUTINE COEFRO
!                    *****************
!
     &(CF,H,U,V,KARMAN,KFROT,CHESTR,GRAV,MESH,T1)
!
!***********************************************************************
! TELEMAC2D   V6P1                                   21/08/2010
!***********************************************************************
!
!brief    COMPUTES THE FRICTION COEFFICIENT CF.
!
!history  J-M HERVOUET (LNHE)
!+        27/07/2009
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
!| CF             |<--| ADIMENSIONAL FRICTION COEFFICIENT
!| CHESTR         |-->| FRICTION COEFFICIENTS
!| GRAV           |-->| GRAVITY
!| H              |-->| WATER DEPTH
!| KARMAN         |-->| VON KARMAN CONSTANT
!| KFROT          |-->| FRICTION LAW ON BOTTOM
!| MESH           |-->| MESH STRUCTURE
!| T1             |<->| WORK BIEF_OBJ STRUCTURE
!| U              |-->| X-COMPONENT OF VELOCITY
!| V              |-->| Y-COMPONENT OF VELOCITY
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN)            :: KFROT
      DOUBLE PRECISION, INTENT(IN)   :: GRAV,KARMAN
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: CF,T1
      TYPE(BIEF_OBJ), INTENT(IN)     :: CHESTR,H,U,V
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER NPOIN,N,IELMC,IELMH
!
      DOUBLE PRECISION TIERS,HC,UNORM,AUX,INLOG,C
      DOUBLE PRECISION, POINTER :: HH(:)
!
      INTRINSIC SQRT,MAX,LOG
!
!-----------------------------------------------------------------------
!
      IELMC = CF%ELM
      IELMH = H%ELM
!
!  DEPTH WITH THE SAME DISCRETISATION AS CF
!  IN CASES WHERE IT IS NEEDED.
!
      IF(KFROT.NE.0.AND.KFROT.NE.2) THEN
!
        IF(IELMC.EQ.IELMH) THEN
          HH=>H%R
        ELSE
          CALL OS( 'X=Y     ' , X=T1 , Y=H )
          CALL CHGDIS( T1 , IELMH , IELMC , MESH )
          HH=T1%R
        ENDIF
!
      ENDIF
!
      NPOIN = CF%DIM1
!
!-----------------------------------------------------------------------
!
      TIERS  = 1.D0/3.D0
!
!  FRICTION COEFFICIENT
!
!     LAWS OF FRICTION:
!
!     KFROT = 0:  NO FRICTION
!     KFROT = 1:  HAALAND
!     KFROT = 2:  CHEZY
!     KFROT = 3:  STRICKLER
!     KFROT = 4:  MANNING
!     KFROT = 5:  NIKURADSE
!
!     *******************
      IF(KFROT.EQ.0) THEN
!     *******************
!
        DO N=1,NPOIN
          CF%R(N) = 0.D0
        ENDDO
!
!     ***********************
      ELSEIF(KFROT.EQ.1) THEN
!     ***********************
!
        DO N=1,NPOIN
          HC = MAX(HH(N),1.D-4)
          UNORM = MAX(SQRT(U%R(N)**2+V%R(N)**2),1.D-6)
!                       1.D-6: LAMINAR VISCOSITY OF THE WATER
          INLOG =(6.9D0*1.D-6/4.D0/HC/UNORM)**3+
     &                                  (CHESTR%R(N)/14.8D0/HC)**3.33
          INLOG = MIN(1.D0-1.D-6,INLOG)
          AUX   = -0.6D0*LOG(INLOG)/LOG(10.D0)
          CF%R(N) = 0.25D0 / AUX**2
        ENDDO
!
!     ***********************
      ELSEIF(KFROT.EQ.2) THEN
!     ***********************
!
        DO N=1,NPOIN
          CF%R(N) = 2 * GRAV / CHESTR%R(N)**2
        ENDDO
!
!     ***********************
      ELSEIF(KFROT.EQ.3) THEN
!     ***********************
!
        DO N=1,NPOIN
          HC = MAX(HH(N),1.D-4)
          CF%R(N) = 2 * GRAV / CHESTR%R(N)**2 / HC**TIERS
        ENDDO
!
!     ***********************
      ELSEIF(KFROT.EQ.4) THEN
!     ***********************
!
        DO N=1,NPOIN
          HC = MAX(HH(N),1.D-4)
          CF%R(N) = 2 * CHESTR%R(N)**2 * GRAV / HC**TIERS
        ENDDO
!
!     ***********************
      ELSEIF(KFROT.EQ.5) THEN
!     ***********************
!
        DO N=1,NPOIN
          HC = MAX(HH(N),1.D-4)/EXP(1.D0)
          CF%R(N) = 2.D0 / (LOG( 30.*HC/CHESTR%R(N))/KARMAN )**2
        ENDDO
!
!     ****
      ELSE
!     ****
!
        IF(LNG.EQ.1) WRITE(LU,300) KFROT
        IF(LNG.EQ.2) WRITE(LU,301) KFROT
300     FORMAT(1X,'COEFRO : LOI DE FROTTEMENT INCONNUE :',1I6)
301     FORMAT(1X,'COEFRO: UNKNOWN LAW OF BOTTOM FRICTION: ',1I6)
        CALL PLANTE(1)
        STOP
!
!     *****
      ENDIF
!     *****
!
!-----------------------------------------------------------------------
!
      RETURN
      END
