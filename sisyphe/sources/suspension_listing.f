      ! ***************************** !
        SUBROUTINE SUSPENSION_LISTING !
      ! ***************************** !

     *(MESH,CST,ZFCL_S,UCONV,VCONV,MASKEL,IELMT,DT,MSK,T1)

C**********************************************************************C
C SISYPHE VERSION 5.8  22/12/04  F. HUVELIN                            C
C**********************************************************************C

                   ! ========================= !
                   ! Listing of min/max values !
                   ! ========================= !


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
C                    <=  Can't be change by the user                   C
C                    =>  Can be changed by the user                    C 
C ---------------------------------------------------------------------C
!                                                                      !
! CALLED BY SUSPENSION_COMPUTATION                                     !
!                                                                      !
! CALL      MAXI   (BIEF)                                              !
!           MINI   (BIEF)                                              !
!           CFLPSI (BIEF)                                              !
!                                                                      !
!======================================================================!
!======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS               !
!======================================================================!
!======================================================================!

      ! 1/ MODULES
      ! ----------
      USE INTERFACE_SISYPHE,
     &    EX_SUSPENSION_LISTING => SUSPENSION_LISTING
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU


      ! 2/ GLOBAL VARIABLES
      ! -------------------
      TYPE(BIEF_MESH),  INTENT(INOUT) :: MESH
      TYPE(BIEF_OBJ),   INTENT(IN)    :: CST,ZFCL_S
      TYPE(BIEF_OBJ),   INTENT(IN)    :: UCONV,VCONV,MASKEL
      INTEGER,          INTENT(IN)    :: IELMT
      DOUBLE PRECISION, INTENT(IN)    :: DT
      LOGICAL,          INTENT(IN)    :: MSK
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: T1


      ! 3/ LOCAL VARIABLES
      ! ------------------
      INTEGER          :: IMAX,IMA
      DOUBLE PRECISION :: XMAX,XMA

      INTEGER                        P_IMAX
      DOUBLE PRECISION P_DMAX,P_DMIN
      EXTERNAL         P_DMAX,P_DMIN,P_IMAX

!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!
!
      CALL MAXI(XMAX,IMAX,CST%R,MESH%NPOIN)
      IF(NCSIZE.GT.1) THEN
        XMA=P_DMAX(XMAX)
        IF(XMAX.EQ.XMA) THEN
          IMA=MESH%KNOLG%I(IMAX)
        ELSE
          IMA=0
        ENDIF
        IMA=P_IMAX(IMA)
      ELSE
        IMA=IMAX
        XMA=XMAX 
      ENDIF
      IF(LNG.EQ.1) WRITE(LU,500) XMA, IMA
      IF(LNG.EQ.2) WRITE(LU,510) XMA, IMA
      CALL MINI(XMAX, IMAX, CST%R, MESH%NPOIN)
      IF(NCSIZE.GT.1) THEN
        XMA=P_DMIN(XMAX)
        IF(XMAX.EQ.XMA) THEN
          IMA=MESH%KNOLG%I(IMAX)
        ELSE
          IMA=0
        ENDIF
        IMA=P_IMAX(IMA)
      ELSE
        IMA=IMAX
        XMA=XMAX 
      ENDIF
      IF(LNG.EQ.1) WRITE(LU,501) XMA, IMA
      IF(LNG.EQ.2) WRITE(LU,511) XMA, IMA
      CALL MAXI(XMAX, IMAX, ZFCL_S%R, MESH%NPOIN)
      IF(NCSIZE.GT.1) THEN
        XMA=P_DMAX(XMAX)
        IF(XMAX.EQ.XMA) THEN
          IMA=MESH%KNOLG%I(IMAX)
        ELSE
          IMA=0
        ENDIF
        IMA=P_IMAX(IMA)
      ELSE
        IMA=IMAX
        XMA=XMAX 
      ENDIF
      IF(LNG.EQ.1) WRITE(LU,502) XMA, IMA
      IF(LNG.EQ.2) WRITE(LU,512) XMA, IMA
      CALL MINI(XMAX, IMAX, ZFCL_S%R, MESH%NPOIN)
      IF(NCSIZE.GT.1) THEN
        XMA=P_DMIN(XMAX)
        IF(XMAX.EQ.XMA) THEN
          IMA=MESH%KNOLG%I(IMAX)
        ELSE
          IMA=0
        ENDIF
        IMA=P_IMAX(IMA)
      ELSE
        IMA=IMAX
        XMA=XMAX 
      ENDIF
      IF(LNG.EQ.1) WRITE(LU,503) XMA, IMA
      IF(LNG.EQ.2) WRITE(LU,513) XMA, IMA
!
      CALL CFLPSI(T1, UCONV, VCONV, DT, IELMT, MESH, MSK, MASKEL)
      CALL MAXI(XMAX, IMAX, T1%R, MESH%NPOIN)
      IF(NCSIZE.GT.1) THEN
        XMA=P_DMAX(XMAX)
        IF(XMAX.EQ.XMA) THEN
          IMA=MESH%KNOLG%I(IMAX)
        ELSE
          IMA=0
        ENDIF
        IMA=P_IMAX(IMA)
      ELSE
        IMA=IMAX
        XMA=XMAX 
      ENDIF
      IF(LNG.EQ.1) WRITE(LU,507) XMA,IMA
      IF(LNG.EQ.2) WRITE(LU,517) XMA,IMA

      !----------------------------------------------------------------!
500   FORMAT(' CONCENTRATION MAXIMALE     : ',G16.7,' %, NOEUD = ',1I8)
501   FORMAT(' CONCENTRATION MINIMALE     : ',G16.7,' %, NOEUD = ',1I8)
502   FORMAT(' EVOLUTION MAXIMALE         : ',G16.7,'  , NOEUD = ',1I8)
503   FORMAT(' EVOLUTION MINIMALE         : ',G16.7,'  , NOEUD = ',1I8)
507   FORMAT(' CFL MAX POUR LA SUSPENSION : ',G16.7,'  , NOEUD = ',1I8)
      !----------------------------------------------------------------!
510   FORMAT(' MAXIMAL CONCENTRATION    : ',G16.7,' %, NODE = ',1I8)
511   FORMAT(' MINIMAL CONCENTRATION    : ',G16.7,' %, NODE = ',1I8)
512   FORMAT(' MAXIMAL EVOLUTION        : ',G16.7,'  , NODE = ',1I8)
513   FORMAT(' MINIMAL EVOLUTION        : ',G16.7,'  , NODE = ',1I8)
517   FORMAT(' MAX. CFL FOR SUSPENSION  : ',G16.7,'  , NODE = ',1I8)
      !----------------------------------------------------------------!

!======================================================================!
!======================================================================!

      RETURN
      END SUBROUTINE SUSPENSION_LISTING
