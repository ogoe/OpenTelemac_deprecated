C                    ******************************** 
                     SUBROUTINE BEDLOAD_SOLIDISCHARGE
C                    ********************************
C
     &(MESH,U2D,V2D,UNORM,HN,TW,UW,MU,TOB,CF,TOBW,FW,THETAW,
     & AVA,MASKPT,MASKEL,ACLADM,UNLADM,KSP,KSR,LIQBOR,
     & QBOR,DEBUG,NPOIN,NPTFR,IELMT,ICF,KENT,OPTBAN,
     & HIDFAC,GRAV,DM,D90,XWC,XMVE,XMVS,XKV,VCE,HMIN,
     & HIDI,KARMAN,ZERO,PI,KARIM_HOLLY_YANG,
     & SUSP,MSK,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,
     & T11,T12,AC,HIDING,QSC,QSS,
     & SLOPEFF,COEFPN,PHISED,CALFA,SALFA,BETA,ZF_C,S, 
     & DEVIA,BETA2,SECCURRENT,BIJK,HOULE,UNSV2D)
C
C**********************************************************************
C SISYPHE VERSION 6.0  15/09/2009  J.-M. HERVOUET                      
C SISYPHE VERSION 5.9  11/03/2008  J.-M. HERVOUET                      
C SISYPHE VERSION 5.5  14/09/2004  F.    HUVELIN                       
C SISYPHE VERSION 5.4  --/10/2003  C.    VILLARET                      
C SISYPHE VERSION 5.3  --/07/2002  M.    GONZALES DE LINARES           
C SISYPHE VERSION 5.2  --/12/2001  B.    MINH DUC                      
C SISYPHE VERSION 5.1  20/05/1995  E.    PELTIER                       
C SISYPHE VERSION 5.1  20/05/1995  C.    LENORMANT                     
C SISYPHE VERSION 5.1  20/05/1995  J.-M. HERVOUET                      
C**********************************************************************
C
C 11/03/2009 : CORRECTIONS POUR LE PARALLELISME 
C
C COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT
C**********************************************************************C
C                                                                      C
C                 SSSS I   SSSS Y   Y PPPP  H   H EEEEE                C
C                S     I  S      Y Y  P   P H   H E                    C
C                 SSS  I   SSS    Y   PPPP  HHHHH EEEE                 C
C                    S I      S   Y   P     H   H E                    C
C                SSSS  I  SSSS    Y   P     H   H EEEEE                C
C                                                                      C
C----------------------------------------------------------------------C
C                             ARGUMENTS                                C
C .________________.____.______________________________________________C
C |      NOM       |MODE|                   ROLE                       C
C |________________|____|______________________________________________C
C |________________|____|______________________________________________C
C                    <=  Can't be changed by the user                  C
C                    =>  Can be changed by the user                    C 
C ---------------------------------------------------------------------C
!                                                                      !
! CALLED BY BEDLOAD_MAIN                                               !
!                                                                      !
! CALL      BEDLOAD_EFFPNT                                             !
!           BEDLOAD_HIDING_FACTOR                                      !
!           BEDLOAD_FORMULA                                            !
!                                                                      !
!                                                                      !
!======================================================================!
!======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS               !
!======================================================================!
!======================================================================!
!
      USE INTERFACE_SISYPHE,
     &    EX_BEDLOAD_SOLIDISCHARGE => BEDLOAD_SOLIDISCHARGE
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C 
      TYPE(BIEF_MESH),  INTENT(INOUT) :: MESH
      TYPE(BIEF_OBJ),   INTENT(IN)    :: U2D, V2D,  HN, TW, UW
      TYPE(BIEF_OBJ),   INTENT(IN)    :: UNORM ,MU, KSR ,KSP  
      TYPE(BIEF_OBJ),   INTENT(IN)    :: TOB, CF, TOBW, FW, THETAW
      TYPE(BIEF_OBJ),   INTENT(IN)    :: MASKPT, MASKEL
      TYPE(BIEF_OBJ),   INTENT(IN)    :: ACLADM, UNLADM, LIQBOR, QBOR
      INTEGER,          INTENT(IN)    :: DEBUG
      INTEGER,          INTENT(IN)    :: NPOIN, NPTFR, IELMT, ICF
      INTEGER,          INTENT(IN)    :: KENT, OPTBAN,HIDFAC
      DOUBLE PRECISION, INTENT(IN)    :: GRAV, DM, D90, XWC, XMVE, XMVS
      DOUBLE PRECISION, INTENT(IN)    :: XKV, VCE, HMIN
      DOUBLE PRECISION, INTENT(IN)    :: HIDI
      DOUBLE PRECISION, INTENT(IN)    :: KARMAN, ZERO, PI
      DOUBLE PRECISION, INTENT(IN)    :: KARIM_HOLLY_YANG
      LOGICAL,          INTENT(IN)    :: SUSP, MSK,SECCURRENT,HOULE
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: T1,T2,T3,T4,T5,T6
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: T7,T8,T9,T10,T11,T12
      DOUBLE PRECISION, INTENT(INOUT) :: AC
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: HIDING
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: QSC,QSS
C
      INTEGER,          INTENT(IN)    :: SLOPEFF,DEVIA 
      DOUBLE PRECISION, INTENT(IN)    :: PHISED,BETA,BETA2 
      TYPE(BIEF_OBJ),   INTENT(IN)    :: ZF_C,S,UNSV2D
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: CALFA,SALFA,COEFPN
C
      DOUBLE PRECISION, INTENT(IN)    :: BIJK,AVA(NPOIN)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION P_DMAX
      EXTERNAL         P_DMAX
C
      INTEGER          :: I
!
!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!
! 
! ******************************************************** 
! 0 - COMPUTATION OF THE PARAMETERS FOR THE SLOPE EFFECT   
! ******************************************************** 
! 
      IF (DEBUG > 0) WRITE(LU,*) 'BEDLOAD_EFFPNT'
C 
C     EFFET DE PENTE 
C 
      IF(DEVIA.EQ.0) THEN
        CALL OS('X=Y/Z   ',CALFA, U2D, UNORM, 0.D0, 2, 1.D0, 1.D-12) 
        CALL OS('X=Y/Z   ',SALFA, V2D, UNORM, 0.D0, 2, 0.D0, 1.D-12) 
      ENDIF
C
      IF(SLOPEFF.EQ.0) CALL OS('X=C     ',X=COEFPN,C=1.D0)
C 
      IF(SLOPEFF.NE.0.OR.DEVIA.NE.0) THEN
C   
      CALL BEDLOAD_EFFPNT 
     &     (MASKEL,LIQBOR,S,ZF_C,U2D,V2D,UNORM,NPOIN,NPTFR,IELMT, 
     &      KENT,BETA,PI,MSK,MESH,T1,T2,T3,T4, 
     &      COEFPN,CALFA,SALFA,SLOPEFF,PHISED,DEVIA,BETA2, 
     &      TOB,XMVS,XMVE,DM,GRAV,UNSV2D) 
C
      ENDIF
C 
      IF (DEBUG > 0) WRITE(LU,*) 'END_BEDLOAD_EFFPNT' 
! 
      ! **************************************** !
      ! I - COEFFICIENT DE MASQUAGE / EXPOSITION ! 
      ! **************************************** !
      IF (DEBUG > 0) WRITE(LU,*) 'BEDLOAD_HIDING_FACTOR'
!
!     WITH HUNZIKER FORMULA (6), THE HIDING FACTOR IS COMPUTED WITH
!     THE SOLID DISCHARGE (SEE BEDLOAD_HUNZ_MEYER.F)
!
      IF(ICF.NE.6) THEN
        CALL BEDLOAD_HIDING_FACTOR
     &     (ACLADM, HIDFAC, NPOIN, HIDI, DM, KARIM_HOLLY_YANG, HIDING)
      ENDIF
      IF (DEBUG > 0) WRITE(LU,*) 'END_BEDLOAD_HIDING_FACTOR'
!
      ! ******************************************* !
      ! II - QSC CALCULATED USING EMPIRICAL FORMULA ! 
      !      T1 = DQSC/DH                           ! 
      ! ******************************************* !
      IF (DEBUG > 0) WRITE(LU,*) 'BEDLOAD_FORMULA'
!
      CALL BEDLOAD_FORMULA
     &  (U2D,V2D, UNORM,HN, CF, MU,TOB, TOBW, UW, TW, THETAW, FW, 
     &   ACLADM, UNLADM, KSP,KSR,AVA, NPOIN, ICF, HIDFAC, XMVS, XMVE,
     &   DM, GRAV, VCE, XKV, HMIN, XWC, D90, KARMAN, ZERO,
     &   PI, SUSP, AC, HIDING, T1, T2, T3, T4, T5, T6, T7, T8, T9, 
     &   T10, T11, T12, QSC, QSS, IELMT,SECCURRENT,
     &   SLOPEFF, COEFPN, BIJK, HOULE)
      IF (DEBUG > 0) WRITE(LU,*) 'END_BEDLOAD_FORMULA'
!
      ! **************************************************** !
      ! IV - TRAITEMENT DES POINTS FRONTIERES A DEBIT IMPOSE ! 
      ! **************************************************** !
      IF (DEBUG > 0) WRITE(LU,*) 'BOUNDARY_NODES_TREATMENT'
      DO I = 1 , NPTFR
        IF(LIQBOR%I(I).EQ.KENT) QSC%R(MESH%NBOR%I(I)) = QBOR%R(I)
      ENDDO
      IF (DEBUG > 0) WRITE(LU,*) 'END_BOUNDARY_NODES_TREATMENT'
!
      ! ************************************ !
      ! V - TRAITEMENT DES BANCS DECOUVRANTS ! 
      ! ************************************ !
      IF(OPTBAN.EQ.2) THEN
        IF (DEBUG > 0) WRITE(LU,*) 'TIDAL_FLATS_TREATMENT'
        CALL OS('X=XY    ', X=QSC, Y=MASKPT)
        IF (DEBUG > 0) WRITE(LU,*) 'END_TIDAL_FLATS_TREATMENT'
      ENDIF
!
!======================================================================!
!======================================================================!
!
      RETURN
      END
