!                    **********************
                     SUBROUTINE TOB_SISYPHE
!                    **********************
!
     & (TOB, TOBW, MU, KS,KSP, KSR,CF,FW,CHESTR,UETCAR,
     &  CF_TEL,KS_TEL,CODE,
     &  KFROT,ICR, KSPRATIO, HOULE,GRAV,XMVE,XMVS, VCE, KARMAN,
     &  ZERO,HMIN,HN, ACLADM, UNORM,UW, TW, NPOIN,KSPRED,IKS)
!
!***********************************************************************
! SISYPHE   V6P1                                   21/07/2011
!***********************************************************************
!
!brief    COMPUTES THE TOTAL STRESS AT THE BOTTOM DEPENDING
!+                ON WHETHER SISYPHE IS COUPLED OR NOT.
!
!history  CV
!+        **/04/05
!+
!+   CORRECTION WHEN SISYPHE IS RUN ALONE: DO NOT MODIFY
!
!history  C. VILLARET (LNHE)
!+        29/11/06
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
!history  C.VILLARET (EDF-LNHE), P.TASSI (EDF-LNHE)
!+        19/07/2011
!+        V6P1
!+  Name of variables   
!+   
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| ACLADM         |-->| MEAN DIAMETER OF SEDIMENT
!| CF             |-->| QUADRATIC FRICTION COEFFICIENT
!| CF_TEL         |-->| QUADRATIC FRICTION COEFFICIENT (COUPLED T2D)
!| CHESTR         |-->| FRICTION COEFFICIENT (KEYWORD)
!| CODE           |-->| CALLING PROGRAM IN COUPLING
!| FW             |-->| QUADRATIC FRICTION COEFFICIENT (WAVE)
!| GRAV           |-->| ACCELERATION OF GRAVITY
!| HMIN           |-->| MINIMUM VALUE OF WATER DEPTH
!| HN             |-->| WATER DEPTH
!| HOULE          |-->| INCLUDE WAVES COMPUTATIONS
!| ICR            |-->| ICR=0: MU=1
!|                |   | ICR=1: SKIN FRICTION CORRECTION USE KSP
!|                |   | ICR=2: RIPPLE ROUGHNESS USE KSR, KSR
!| KARMAN         |-->| VON KARMAN CONSTANT  
!| KFROT          |-->| FRICTION LAW
!| KS             |<--| RUGOSITE TOTALE
!| KSP            |<--| RUGOSITE DE PEAU
!| KSPRATIO       |-->| RATIO BETWEEN SKIN BED ROUGHNESS AND GRAIN DIAMETER
!| KSR            |<--| RUGOSITE DE RIDE
!| KS_TEL         |<--| RUGOSITE TOTALE
!| MU             |<->| CORRECTION FACTOR FOR BED ROUGHNESS
!| NPOIN          |-->| NUMBER OF POINTS
!| TOB            |<->| BED SHEAR STRESS (TOTAL FRICTION)
!| TOBW           |-->| WAVE INDUCED SHEAR STRESS
!| TW,UW          |-->| WAVE PERIOD AND ORBITAL VELOCITY
!| UETCAR         |-->| SQUARE OF THE FRICTION VELOCITY (COUPLED T3D)
!| UNORM          |-->| INTENSITE DU COURANT
!| VCE            |-->| WATER VISCOSITY
!| XMVE           |-->| FLUID DENSITY (MASS) 
!| XMVS           |-->| SEDIMENT DENSITY (MASS) 
!| ZERO           |-->| ZERO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
      USE INTERFACE_SISYPHE, EX_TOB_SISYPHE=>TOB_SISYPHE
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER,            INTENT(IN)  :: NPOIN,KFROT,ICR, IKS
      LOGICAL,            INTENT(IN)  :: KSPRED
      LOGICAL,            INTENT(IN)  :: HOULE
      CHARACTER(LEN=24),  INTENT(IN)  :: CODE
      DOUBLE PRECISION,   INTENT(IN)  :: XMVE,XMVS, VCE,GRAV,KARMAN
      DOUBLE PRECISION,   INTENT(IN)  :: ZERO,HMIN,KSPRATIO
      TYPE(BIEF_OBJ), INTENT(IN)      :: UETCAR
      TYPE(BIEF_OBJ), INTENT(IN)      :: HN,UNORM
      TYPE(BIEF_OBJ), INTENT(IN)      :: TW,UW
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: KS,KSP,KSR,KS_TEL
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: CHESTR,MU
      TYPE(BIEF_OBJ), INTENT(IN)      :: ACLADM
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: CF,TOB
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: FW,TOBW
      TYPE(BIEF_OBJ), INTENT(IN)      :: CF_TEL
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER                     :: I
      DOUBLE PRECISION            :: A,B,C, HCLIP,KSMAX
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! ----------------------------------------------------------------------------------------------
!  QUADRATIC FRICTION COEFICIENT       :  ---> CF
!-----------------------------------------------------------------------
!
!     INTERNAL COUPLING WITH TELEMAC2D OR 3D
!         	UETCAR IS CF IN TELEMAC-2D
!         	UETCAR IS UETCAR IN TELEMAC3D ?
!  KSP : skin friction
!  KSR: ripple roughness
!  KS : total bed roughness
!  initialisation
!
      CALL OS('X=CY    ', X=KSP, Y=ACLADM, C=KSPRATIO)
      CALL OS('X=CY    ', X=KSR, Y=ACLADM, C=KSPRATIO)
!
      IF(KSPRED) THEN
!
!        bed roughness predictor
!
         CALL KS_SISYPHE(IKS,KS,KSP,KSR,KSPRATIO,HOULE,
     &                   GRAV,XMVE,XMVS,VCE,
     &                   HMIN,HN,ACLADM,UNORM,UW,TW,NPOIN)
         CALL COEFRO_SISYPHE(CF,HN,KFROT,KS,GRAV,NPOIN,HMIN,KARMAN)
         IF(CODE(1:7).EQ.'TELEMAC') CALL OS('X=Y     ', X=KS_TEL, Y=KS)
!
      ELSE
!
! here the total bed roughness is calculated as a function of friction coefficient
! -- > issued from Telemac if coupling
! -- > from the steering file of Sisyphe
!  
        IF(CODE(1:7).EQ.'TELEMAC') THEN
          CALL OV('X=Y     ',CF%R,CF_TEL%R,CF_TEL%R,0.D0,CF%DIM1)
        ELSE
          CALL COEFRO_SISYPHE(CF,HN,KFROT,CHESTR,GRAV,NPOIN,HMIN,KARMAN)
        ENDIF
        DO I =1,NPOIN
          A = -KARMAN*SQRT(2.D0/MAX(CF%R(I),ZERO))
          KS%R(I)=12.D0*HN%R(I)*EXP(A)
          KS%R(I)=MAX(KS%R(I),KSP%R(I))
        ENDDO
!
      ENDIF
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Frottement total: loi quadratique sauf couplage 3D
!  --> TOB
!
!     INTERNAL COUPLING WITH TELEMAC3D
!     UETCAR CORRESPONDS TO THE FRICTION VELOCITY SQUARED
!
      IF(CODE(1:9).EQ.'TELEMAC3D') THEN
        CALL OS( 'X=CY     ',X=TOB,Y=UETCAR,C=XMVE)
      ELSE
        DO I=1,NPOIN
          TOB%R(I) = XMVE*0.5D0*CF%R(I)*UNORM%R(I)**2
        ENDDO
      ENDIF
!
! -----WAVE-INDUCED FRICTION -----------------------------
!  --> TOBW
!
      IF(HOULE) THEN
        CALL TOBW_SISYPHE(TOBW%R,CF%R,FW%R,UW%R,TW%R,HN%R,NPOIN,XMVE)
      ENDIF
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! SKIN FRICTION CORRECTOR
!                ---> MU = TOP/TOB
! ICR=0:    MU=1
! ICR=1     : SKIN FRICTION CORRECTION USE KSP
! ICR= 2    : RIPPLE ROUGHNESS USE KSR, KSR
! COUPLED WITH TELEMAC: MU>1 IS ACCEPTABLE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      IF(ICR.EQ.0) THEN
        CALL OS('X=C     ', X=MU, C=1.D0)
      ELSEIF(ICR.EQ.1) THEN
        DO I= 1, NPOIN
          IF(CF%R(I).GT.ZERO.AND.HN%R(I).GT.KSP%R(I)) THEN
            HCLIP=MAX(HN%R(I),KSP%R(I))
            A = 2.5D0*LOG(12.D0*HCLIP/KSP%R(I))
            C =2.D0/A**2
            MU%R(I) = C/CF%R(I)
          ELSE
            MU%R(I) = 0.D0
          ENDIF
        ENDDO
      ELSEIF(ICR.EQ.2) THEN
        DO I= 1, NPOIN
          KSMAX=MAX(KSR%R(I),KSP%R(I))
          IF(HN%R(I).GT.KSMAX.AND.CF%R(I).GT.ZERO)THEN
            HCLIP=MAX(HN%R(I),KSMAX)
            A = LOG(12.D0*HCLIP/KSP%R(I))
            B = LOG(12.D0*HCLIP/KSR%R(I))
            C = 0.32D0/CF%R(I)
            MU%R(I) = C/SQRT(B*A**3)
          ELSE
            MU%R(I) = 0.D0
          ENDIF
        ENDDO
      ENDIF
!
!------------------------------------------------------------
!
      RETURN
      END
