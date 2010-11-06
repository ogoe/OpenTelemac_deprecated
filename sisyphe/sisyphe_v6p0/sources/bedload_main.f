      ! *********************** !
        SUBROUTINE BEDLOAD_MAIN !
      ! *********************** !

     &  (ACLADM,KSP,KSR, V2DPAR,UNSV2D,CF,EBOR,FW,HN,LIQBOR, 
     &   MASK, MASKEL, MASKPT, Q, QBOR, U2D,
     &   V2D, S,UNLADM,UW,THETAW,MU,TOB,TOBW,TW,ZF,
     &   DEBUG, HIDFAC, ICF, IELMT, ISOUS, KDDL, KDIR,
     &   KENT, KINC, KLOG, KNEU, KSORT, LOADMETH, LT,
     &   NPOIN, NPTFR, NSICLA, OPTBAN, LS0, BETA, FD90, FDM,
     &   GRAV, HIDI, HMIN, VCE, XKV, XMVE, XMVS, XWC,
     &   PI, KARMAN, ZERO, KARIM_HOLLY_YANG,MSK, SUSP, VF,
     &   ENTET, CONST_ALAYER, LCONDIS, LGRAFED, MESH,
     &   ELAY, LIEBOR, LIMTEC, MASKTR,
     &   IT1, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11,
     &   T12,T13,UNORM,AC, AT0, DTS, ELAY0, FRACSED_GF,
     &   AVAIL, BREACH, CALFA, COEFPN,
     &   DZF_GF, HIDING, QSCL_C, QSCL_S, QS_C,
     &   QSCLXC, QSXC, QSCLYC, QSYC, SALFA, ZF_C, ZFCL_C, NSOUS,
     &   ENTETS, SECCURRENT, SLOPEFF, 
     &   PHISED, DEVIA, BETA2, BIJK,SEDCO,HOULE)

C**********************************************************************
C SISYPHE VERSION 6.0                              14/09/04  F. HUVELIN
C
C JMH LE 21/12/2006: BEDLOAD_TIMESTEP SUPPRIME
C**********************************************************************


             ! ========================================== !
             ! Main subroutine for the bed-load transport !
             ! ========================================== !


C COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT
C**********************************************************************
C                                                                      
C                 SSSS I   SSSS Y   Y PPPP  H   H EEEEE                
C                S     I  S      Y Y  P   P H   H E                    
C                 SSS  I   SSS    Y   PPPP  HHHHH EEEE                 
C                    S I      S   Y   P     H   H E                    
C                SSSS  I  SSSS    Y   P     H   H EEEEE                
C                                                                      
C----------------------------------------------------------------------C
C                             ARGUMENTS                                C
C .________________.____.______________________________________________C
C |      NOM       |MODE|                   ROLE                       C
C |________________|____|______________________________________________C
C |________________|____|______________________________________________C
C                    <=  Can't be change by the user                   C
C                    =>  Can be changed by the user                    C 
C ---------------------------------------------------------------------C
!                                                                      !
! CALLED BY SISYPHE                                                    !
!                                                                      !
! CALL      BEDLOAD_INIT                                               !
!           BEDLOAD_SOLIDISCHARGE                                      !
!           BEDLOAD_EVOL                                               !
!           BEDLOAD_POSTREATMENT                                       !
!                                                                      !
!======================================================================!
!======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS               !
!======================================================================!
!======================================================================!
!
      ! 1/ MODULES
      ! ----------
      USE BIEF
      USE INTERFACE_SISYPHE, EX_BEDLOAD_MAIN => BEDLOAD_MAIN
      USE DECLARATIONS_SISYPHE, ONLY : DREDGESIM
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU


      ! 2/ GLOBAL VARIABLES
      ! -------------------

      TYPE(BIEF_OBJ),   INTENT(IN)    :: ACLADM, KSR,V2DPAR,UNSV2D
      TYPE(BIEF_OBJ),   INTENT(IN)    :: CF,FW,KSP,HN,LIQBOR
      TYPE(BIEF_OBJ),   INTENT(IN)    :: MASK, MASKEL, MASKPT
      TYPE(BIEF_OBJ),   INTENT(IN)    :: Q, QBOR
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: U2D, V2D,TOB, MU,UNORM,EBOR
      TYPE(BIEF_OBJ),   INTENT(IN)    :: S,UNLADM
      TYPE(BIEF_OBJ),   INTENT(IN)    :: UW, THETAW,  TOBW, TW
      TYPE(BIEF_OBJ),   INTENT(IN)    :: ZF
      INTEGER,          INTENT(IN)    :: DEBUG, HIDFAC, ICF
      INTEGER,          INTENT(IN)    :: IELMT, ISOUS, KDDL, KDIR, KENT
      INTEGER,          INTENT(IN)    :: KINC, KLOG, KNEU, KSORT
      INTEGER,          INTENT(IN)    :: LOADMETH, LT,NPOIN, NPTFR
      INTEGER,          INTENT(IN)    :: NSICLA, OPTBAN
      DOUBLE PRECISION, INTENT(IN)    :: LS0, BETA, FD90(NSICLA)
      DOUBLE PRECISION, INTENT(IN)    :: FDM(NSICLA),GRAV
      DOUBLE PRECISION, INTENT(IN)    :: HIDI(NSICLA),HMIN,VCE
      DOUBLE PRECISION, INTENT(IN)    :: XKV,XMVE,XMVS,XWC(NSICLA)
      DOUBLE PRECISION, INTENT(IN)    :: PI,KARMAN,ZERO
      DOUBLE PRECISION, INTENT(IN)    :: KARIM_HOLLY_YANG
      LOGICAL,          INTENT(IN)    :: MSK, SUSP, VF
      LOGICAL,          INTENT(IN)    :: ENTET, CONST_ALAYER
      LOGICAL,          INTENT(IN)    :: LCONDIS, LGRAFED,SECCURRENT
      LOGICAL,          INTENT(IN)    :: SEDCO(NSICLA),HOULE
      TYPE(BIEF_MESH),  INTENT(INOUT) :: MESH
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: ELAY
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: LIEBOR, LIMTEC, MASKTR
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: IT1,T1,T2,T3,T4,T5,T6,T7
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: T8,T9,T10,T11,T12,T13
      DOUBLE PRECISION, INTENT(INOUT) :: AC(NSICLA), AT0, DTS, ELAY0
      DOUBLE PRECISION, INTENT(INOUT) :: FRACSED_GF(NSICLA)
      DOUBLE PRECISION, INTENT(INOUT) :: AVAIL(NPOIN,10,NSICLA)
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: BREACH, CALFA, COEFPN
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: DZF_GF
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: HIDING
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: QSCL_C,QSCL_S
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: QS_C, QSCLXC, QSXC, QSCLYC
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: QSYC, SALFA, ZF_C, ZFCL_C
      INTEGER,          INTENT(INOUT) :: NSOUS
      LOGICAL,          INTENT(INOUT) :: ENTETS
      DOUBLE PRECISION,   INTENT(IN)  :: BETA2, PHISED
      INTEGER, INTENT (IN)            :: SLOPEFF, DEVIA 
      DOUBLE PRECISION, INTENT(IN)    :: BIJK
  
      ! 3/ LOCAL VARIABLES
      ! ------------------
      INTEGER :: I
!
!
!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!
! 
! 
!     INITIALIZING TECHNICAL BOUNDARIES CONDITIONS
!     
      IF (DEBUG > 0) WRITE(LU,*) 'BEDLOAD_DIFFIN' 
      CALL BEDLOAD_DIFFIN 
     &        (U2D, V2D, MESH%NBOR, MESH%XNEBOR, MESH%YNEBOR, 
     &         MESH%KP1BOR, 
     &         MASKEL, MESH%NELBOR, NPTFR, KENT, KSORT, KLOG, KINC, 
     &         KDIR, KDDL, KNEU, MSK, IT1, LIEBOR, MASKTR, LIMTEC) 
      IF (DEBUG > 0) WRITE(LU,*) 'END_BEDLOAD_DIFFIN' 
!
      DO I = 1, NSICLA
!
!        FOR SAND
         IF(.NOT.SEDCO(I)) THEN
           IF (DEBUG > 0) WRITE(LU,*)
     &       'BEDLOAD_SOLIDISCHARGE : ',I,'/',NSICLA
!   
           CALL BEDLOAD_SOLIDISCHARGE
     &        (MESH, U2D, V2D, UNORM,HN, TW, UW, MU,TOB,CF, 
     &          TOBW,FW,THETAW,AVAIL(1:NPOIN,1,I),
     &          MASKPT, MASKEL, ACLADM,
     &          UNLADM,KSP,KSR, LIQBOR, QBOR%ADR(I)%P, DEBUG, NPOIN,
     &          NPTFR, IELMT, ICF, KENT, OPTBAN, HIDFAC, GRAV,
     &          FDM(I), FD90(I), XWC(I), XMVE, XMVS, XKV, VCE, HMIN,
     &          HIDI(I),KARMAN,ZERO,PI,
     &          KARIM_HOLLY_YANG,SUSP,MSK,T1,T2,
     &          T3, T4, T5, T6, T7, T8, T9, T10, T11,T12, AC(I),
     &          HIDING,QSCL_C%ADR(I)%P,QSCL_S%ADR(I)%P,
     &          SLOPEFF,COEFPN,PHISED,CALFA,SALFA,BETA,ZF,S,  
     &          DEVIA, BETA2 , SECCURRENT, BIJK,HOULE,UNSV2D)
          IF(DEBUG > 0) WRITE(LU,*) 'END_BEDLOAD_SOLIDISCHARGE'
        ELSE
!         FOR COHESIVE SEDIMENT: ZERO bedload transport rate 
!         JMH: IS THIS USEFUL ???         
          CALL OS('X=0     ',X=QSCL_C%ADR(I)%P)
          CALL OS('X=0     ',X=QSCLXC%ADR(I)%P)
          CALL OS('X=0     ',X=QSCLYC%ADR(I)%P)
        ENDIF
C
      ENDDO
!
!     EVOLUTION COMPUTATION FOR EACH CLASS
!
      DO I = 1, NSICLA
!
        IF(.NOT.SEDCO(I)) THEN       
!
          IF (DEBUG > 0) WRITE(LU,*) 'BEDLOAD_EVOL : ',I,'/',NSICLA
          CALL BEDLOAD_EVOL(HN,Q,S,ELAY,ACLADM,AVAIL(1:NPOIN,1,I),
     &                      COEFPN,CALFA,SALFA,LIMTEC,
     &                      EBOR%ADR(I)%P,MASKEL,MASK,
     &                      V2DPAR,UNSV2D,DEBUG,NPOIN,NPTFR,IELMT,
     &                      KENT,KDIR,KDDL,LOADMETH,
     &                      DTS,FDM(I),FD90(I),HMIN,LS0,GRAV,XMVS,XMVE,
     &                      VCE,VF,ENTETS,MSK,CONST_ALAYER,
     &                      LCONDIS,MESH,QSCL_C%ADR(I)%P,
     &                      T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,
     &                      T13,ELAY0,BREACH,QSCLXC%ADR(I)%P,
     &                      QSCLYC%ADR(I)%P,ZFCL_C%ADR(I)%P,SLOPEFF)
          IF(DEBUG.GT.0) WRITE(LU,*) 'END_BEDLOAD_EVOL'
C
        ELSE
C
C         NO EVOLUTION FOR COHESIVE SEDIMENT
          CALL OS('X=0     ',X=ZFCL_C%ADR(I)%P)             
C
        ENDIF
C
      ENDDO
C
C     CALLING DREDGESIM
C
      IF(DREDGESIM) CALL DREDGESIM_INTERFACE(2)

      ! *********************************************** !
      ! II - EVOLUTIONS AND QS FOR EACH CLASS ARE ADDED ! 
      ! *********************************************** !

      ! II.1 - INITIALIZATION 
      ! ---------------------
      
      CALL OS('X=0     ', X=QS_C)
      CALL OS('X=0     ', X=ZF_C)

      ! II.2 - ADD THE CLASSES
      ! ----------------------
      DO I=1,NSICLA
         ! Correction du transport solide pour ne pas prendre en compte
         ! le coef. de porosite (il doit etre non nul...)
         ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IF(.NOT.SEDCO(I)) THEN
          CALL OS('X=CX    ', X=QSCL_C%ADR(I)%P, C=1.D0/XKV)
          CALL OS('X=CX    ', X=QSCLXC%ADR(I)%P, C=1.D0/XKV)
          CALL OS('X=CX    ', X=QSCLYC%ADR(I)%P, C=1.D0/XKV)
          CALL OS('X=X+Y   ', X=QS_C, Y=QSCL_C%ADR(I)%P)
          CALL OS('X=X+Y   ', X=ZF_C, Y=ZFCL_C%ADR(I)%P)
        ENDIF
      ENDDO          
!
!     TIDAL FLATS WITH MASKING     JMH ON 27/07/2006
!
      IF(OPTBAN.EQ.2) CALL OS('X=XY    ',X=ZF_C,Y=MASKPT)
!
      ! II.3 - SLOPE EFFECT FOR THE SUM OF THE QS 
      ! -----------------------------------------
      ! QS : COEFPN already added in QSCL_C
      CALL OS('X=YZ    ', X=QSXC, Y=QS_C, Z=CALFA)
      CALL OS('X=YZ    ', X=QSYC, Y=QS_C, Z=SALFA)

!======================================================================!
!======================================================================!

      RETURN
      END
