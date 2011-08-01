      ! *********************************** !
        MODULE INTERFACE_SISYPHE_SUSPENSION !
        INTERFACE                           !
      ! *********************************** !

C**********************************************************************C
C SISYPHE VERSION 5.9  25/06/2008                             F. HUVELIN
C**********************************************************************C

             ! ======================================= !
             ! Interface for the subroutine of Sisyphe !
             !       for the suspension transport      !
             ! ======================================= !


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

!======================================================================!
!======================================================================!


      !================================================================!
      SUBROUTINE SUSPENSION_BILAN 
     &(MESH,CST,HN,ZFCL_S,MASKEL,IELMT,ITRA,LT,NIT,DT,CSF,
     & MASSOU,MASED0,MSK,ENTET,MASTEN,MASTOU,MASINI,T2,
     & T3,MASFIN,MASDEPT,MASDEP,AGGLOT,
     & VOLU2D,NUMLIQ,NFRLIQ,NPTFR,FLBORTRA)
      USE BIEF_DEF
      IMPLICIT NONE
      TYPE(BIEF_MESH),  INTENT(INOUT) :: MESH
      TYPE(BIEF_OBJ),   INTENT(IN)    :: CST,HN,VOLU2D
      TYPE(BIEF_OBJ),   INTENT(IN)    :: ZFCL_S,MASKEL,FLBORTRA
      INTEGER,          INTENT(IN)    :: IELMT,ITRA,LT,NIT
      INTEGER,          INTENT(IN)    :: NFRLIQ,NPTFR
      INTEGER,          INTENT(IN)    :: NUMLIQ(NFRLIQ)
      DOUBLE PRECISION, INTENT(IN)    :: DT,CSF
      DOUBLE PRECISION, INTENT(IN)    :: MASSOU,MASED0,AGGLOT
      LOGICAL,          INTENT(IN)    :: MSK,ENTET
      DOUBLE PRECISION, INTENT(INOUT) :: MASTEN,MASTOU,MASINI
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: T2,T3
      DOUBLE PRECISION, INTENT(INOUT) :: MASFIN,MASDEPT,MASDEP
      END SUBROUTINE SUSPENSION_BILAN
      !================================================================!

      !================================================================!
      SUBROUTINE SUSPENSION_COMPUTATION
      ! ********************************* !
     &  (SLVTRA, HN,HN_TEL,UCONV, VCONV,  CF, MU,TOB,ACLADM, KSP,KSR,
     &   ELAY, AVA, AFBOR, BFBOR, LIMDIF, CLT, MASKEL, MASKTR,
     &   MASKPT, IFAMAS, NPOIN, IELMT, NPTFR, ITRA, LT, NIT, RESOL,
     &   OPTBAN, KENT,KDDL,KDIR,KSORT,KLOG,KINC,KNEU,
     &   OPTSUP, OPDTRA, DEBUG, CSF_VASE, CSF_SABLE,
     &   TETA_SUSP, DT, MASED0, ZERO, XWC, KARMAN, XMVE, XMVS, GRAV,
     &   HMIN, VITCD, VITCE,PARTHENIADES, ENTETS,
     &   ENTET, BILMA, MSK, CHARR, IMP_INFLOW_C, MESH, ZF, CS,
     &   CST,CTILD,CBOR,DISP,IT1,IT2,IT3,IT4,TB,T1,T2,T3,
     &   T4, T5, T6, T7, T8, T9, T10, T11, T12, W1, TE1, TE2, TE3, S,
     &   AM1_S, AM2_S, MBOR,MASTEN, MASTOU, MASINI, AC,
     &   ZFCL_S, FLUDPT, FLUDP, FLUER, HPROP, DISP_C, CSTAEQ,
     &   MASFIN, MASDEPT, MASDEP, MASSOU,QS_C,ICQ,ZREF,
     &   CORR_CONV,U2D,V2D,SEDCO,DIFT,
     &   DM1,ZCONV,UCONV_TEL,VCONV_TEL,SOLSYS,FLBOR_TEL,FLBOR_SIS,
     &   FLBORTRA,CODE,
     &   VOLU2D,V2DPAR,UNSV2D,NUMLIQ,NFRLIQ,LICBOR,MIXTE,AVAIL,NSICLA,
     &   ES,NCOUCH_TASS,CONC_VASE,TOCE_VASE,
     &   FLUER_VASE,TOCE_MIXTE,MS_SABLE,MS_VASE,TASS)
      USE BIEF_DEF
      IMPLICIT NONE
      TYPE (SLVCFG),    INTENT(INOUT) :: SLVTRA
      TYPE (BIEF_OBJ),  INTENT(IN)    :: ZF,VOLU2D,V2DPAR,UNSV2D
      TYPE (BIEF_OBJ),  INTENT(IN), TARGET    :: HN,HN_TEL
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: UCONV, VCONV
      TYPE (BIEF_OBJ),  INTENT(IN)    :: CF,MU,TOB,ACLADM,KSP,KSR
      TYPE (BIEF_OBJ),  INTENT(IN)    :: AFBOR,BFBOR,ELAY,LICBOR
      TYPE (BIEF_OBJ),  INTENT(IN)    :: MASKEL,MASKPT,IFAMAS
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: MASKTR,LIMDIF,CLT
      INTEGER,          INTENT(IN)    :: NPOIN, IELMT, NPTFR, ITRA, LT
      INTEGER,          INTENT(IN)    :: NIT, RESOL,OPTBAN, KENT,KDDL
      INTEGER,          INTENT(IN)    :: KDIR, OPTSUP, OPDTRA,SOLSYS
      INTEGER,          INTENT(IN)    :: KSORT,KLOG,KINC,KNEU
      INTEGER,          INTENT(IN)    :: NFRLIQ,NSICLA,NCOUCH_TASS
      INTEGER,          INTENT(IN)    :: DEBUG
      INTEGER,          INTENT(IN)    :: NUMLIQ(NFRLIQ)
      DOUBLE PRECISION, INTENT(IN)    :: CSF_VASE, TETA_SUSP, DT, MASED0
      DOUBLE PRECISION, INTENT(IN)    :: ZERO, XWC,CSF_SABLE,AVA(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: KARMAN, XMVE, XMVS, GRAV, HMIN
      DOUBLE PRECISION, INTENT(IN)    :: VITCD, VITCE,PARTHENIADES
      LOGICAL,          INTENT(IN)    :: ENTETS,ENTET,BILMA,MSK,SEDCO
      LOGICAL,          INTENT(IN)    :: CHARR, IMP_INFLOW_C,CORR_CONV
      LOGICAL,          INTENT(IN)    :: DIFT,MIXTE,TASS
      TYPE (BIEF_MESH), INTENT(INOUT) :: MESH
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: CS,CST,CTILD,CBOR,FLBOR_SIS
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: DISP,IT1,IT2,IT3,IT4,TB
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: T2, T3, T4, T5, T6, T7, T8
      TYPE (BIEF_OBJ),  INTENT(INOUT), TARGET :: T1
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: T9, T10, T11, T12, W1, TE1
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: TE2, TE3, S, AM1_S, AM2_S
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: MBOR,ZREF
      DOUBLE PRECISION, INTENT(INOUT) :: MASTEN, MASTOU, MASINI, AC
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: ZFCL_S
      TYPE (BIEF_OBJ),  INTENT(IN)    :: UCONV_TEL,VCONV_TEL
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: FLUDPT,FLUDP,FLUER,FLBORTRA
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: HPROP, DISP_C, CSTAEQ
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: FLUER_VASE,TOCE_MIXTE
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: MS_SABLE,MS_VASE
      DOUBLE PRECISION, INTENT(INOUT) :: MASFIN, MASDEPT, MASDEP
      DOUBLE PRECISION, INTENT(INOUT) :: MASSOU,AVAIL(NPOIN,10,NSICLA)
      DOUBLE PRECISION, INTENT(INOUT) :: ES(NPOIN,10)
      DOUBLE PRECISION, INTENT(INOUT) :: CONC_VASE(10),TOCE_VASE(10)
      TYPE (BIEF_OBJ),  INTENT(IN)    :: QS_C,U2D,V2D,DM1,ZCONV
      TYPE (BIEF_OBJ),  INTENT(IN)    :: FLBOR_TEL
      INTEGER,          INTENT(IN)    :: ICQ
      CHARACTER(LEN=24), INTENT(IN)   :: CODE      
      !----------------------------------------------------------------!
      END SUBROUTINE SUSPENSION_COMPUTATION
      !================================================================!

      !================================================================!
        SUBROUTINE SUSPENSION_CONV ! 
      !================================================================!
     * (TOB,XMVE,CF,NPOIN,ZREF,U2D,V2D,HN,HMIN,
     *     UCONV,VCONV,KARMAN,ZERO,XWC,T1,ALPHA)
      USE BIEF_DEF
      IMPLICIT NONE         
      TYPE (BIEF_OBJ),  INTENT(IN)    :: HN,U2D,V2D,ZREF
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: UCONV, VCONV,T1,ALPHA
      TYPE (BIEF_OBJ),  INTENT(IN)    :: CF, TOB
      INTEGER,          INTENT(IN)    :: NPOIN
      DOUBLE PRECISION, INTENT(IN)    :: ZERO, XWC,HMIN
      DOUBLE PRECISION, INTENT(IN)    :: KARMAN, XMVE          
      !----------------------------------------------------------------!
      END SUBROUTINE SUSPENSION_CONV
      !================================================================!

      !================================================================!
      SUBROUTINE SUSPENSION_DISPERSION
      !----------------------------------------------------------------!
     & (TOB, XMVE,HN,  OPTDIF, NPOIN, XKX, XKY,
     &   T1, T2, T3, KX, KY, KZ, DISP,U2D,V2D,VISC_TEL,CODE)
      !----------------------------------------------------------------!
      USE BIEF
      TYPE (BIEF_OBJ),  INTENT(IN)    :: TOB, HN,VISC_TEL
      INTEGER,          INTENT(IN)    :: OPTDIF, NPOIN
      DOUBLE PRECISION, INTENT(IN)    :: XMVE, XKX, XKY
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: T1, T2, T3
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: KX, KY, KZ, DISP
      TYPE (BIEF_OBJ),  INTENT(IN)    :: U2D,V2D
      CHARACTER(LEN=24), INTENT(IN)   :: CODE
      !----------------------------------------------------------------!
      END SUBROUTINE SUSPENSION_DISPERSION
      !================================================================!

      !==============================!
        SUBROUTINE SUSPENSION_DEPOT  !
      ! ***************************  !

     &(TOB, HN,ACLADM,NPOIN, HMIN,XWC,VITCD,
     & ZERO,KARMAN,XMVE,T1,T2,ZREF,FLUDPT,DEBUG,SEDCO)
      USE BIEF_DEF
      IMPLICIT NONE
      TYPE (BIEF_OBJ),  INTENT(IN)    ::  TOB,HN, ACLADM
      INTEGER,          INTENT(IN)    ::  NPOIN,DEBUG
      LOGICAL,          INTENT(IN)    :: SEDCO
      DOUBLE PRECISION, INTENT(IN)    ::  HMIN
      DOUBLE PRECISION, INTENT(IN)    :: XWC
      DOUBLE PRECISION, INTENT(IN)    :: VITCD
      DOUBLE PRECISION, INTENT(IN)    :: ZERO, KARMAN,XMVE
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: T1, T2
      TYPE (BIEF_OBJ),  INTENT(IN)    :: ZREF
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: FLUDPT
      !----------------------------------------------------------------!
      END SUBROUTINE SUSPENSION_DEPOT
      !================================================================!

      !================================!
        SUBROUTINE SUSPENSION_EROSION  !
      ! ****************************** !

     &(TAUP,HN,ACLADM,AVA,NPOIN,CHARR,XMVE,XMVS,GRAV,HMIN,XWC,
     & ZERO,ZREF,AC,FLUER,CSTAEQ,QSC,ICQ,DEBUG)
!
      USE BIEF_DEF
      IMPLICIT NONE
      TYPE (BIEF_OBJ),  INTENT(IN)    :: TAUP,HN,ACLADM
      INTEGER,          INTENT(IN)    :: NPOIN,DEBUG
      LOGICAL,          INTENT(IN)    :: CHARR
      DOUBLE PRECISION, INTENT(IN)    :: XMVE,XMVS,GRAV,HMIN
      DOUBLE PRECISION, INTENT(IN)    :: XWC,ZERO,AVA(NPOIN)
      TYPE (BIEF_OBJ),  INTENT(IN)    :: ZREF
      DOUBLE PRECISION, INTENT(INOUT) :: AC
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: FLUER,CSTAEQ
      TYPE(BIEF_OBJ),   INTENT(IN)    :: QSC
      INTEGER,          INTENT (IN)   :: ICQ
      !----------------------------------------------------------------!
      END SUBROUTINE SUSPENSION_EROSION
      !================================================================!

      !================================!
       SUBROUTINE SUSPENSION_EROSION_COH  
      ! ****************************** !
     &(TAUP,NPOIN,
     & XMVE,XMVS,GRAV, VITCE,
     & PARTHENIADES,ZERO,DEBUG, 
     & FLUER, ES, TOCE_VASE, NCOUCH_TASS, DT, MS_VASE,TASS)
!
      USE BIEF_DEF
      IMPLICIT NONE 
      TYPE (BIEF_OBJ),  INTENT(IN)    :: TAUP
      INTEGER,          INTENT(IN)    :: NPOIN,DEBUG
      DOUBLE PRECISION, INTENT(IN)    :: XMVE,XMVS,GRAV
      DOUBLE PRECISION, INTENT(IN)    :: VITCE
      DOUBLE PRECISION, INTENT(IN)    :: ZERO,PARTHENIADES 
      DOUBLE PRECISION,  INTENT(INOUT) :: MS_VASE(NPOIN,10)
      DOUBLE PRECISION, INTENT(IN)     :: TOCE_VASE(10), DT
      INTEGER,          INTENT(IN)    :: NCOUCH_TASS
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: FLUER
      DOUBLE PRECISION, INTENT(INOUT) :: ES(NPOIN,10)
      LOGICAL, INTENT(IN) :: TASS
!      
      !----------------------------------------------------------------!
      END SUBROUTINE SUSPENSION_EROSION_COH
      !================================================================!

      !================================!
        SUBROUTINE SUSPENSION_FLUX_MIXTE 
      ! ****************************** !
     &  (TAUP,HN,ACLADM,CS,NPOIN,
     &   CHARR,XMVE,XMVS,GRAV,HMIN,XWC,
     &   ZERO,KARMAN,PARTHENIADES,FLUER_SABLE,FLUER_VASE,ZREF,
     &   AC,CSTAEQ,QSC,ICQ,DEBUG,AVAIL,NSICLA,ES,
     &   TOCE_VASE,NCOUCH_TASS,DT,TOCE_MIXTE,MS_SABLE,MS_VASE)

      USE BIEF_DEF
      IMPLICIT NONE
      TYPE (BIEF_OBJ),  INTENT(IN)    :: TAUP,HN,ACLADM,CS
      INTEGER,          INTENT(IN)    :: NPOIN,DEBUG,NSICLA
      INTEGER,          INTENT(IN)    :: NCOUCH_TASS
      LOGICAL,          INTENT(IN)    :: CHARR
      DOUBLE PRECISION, INTENT(IN)    :: XMVE, XMVS, GRAV, HMIN
      DOUBLE PRECISION, INTENT(IN)    ::  XWC
      DOUBLE PRECISION, INTENT(IN)    :: ZERO, KARMAN, PARTHENIADES
      TYPE (BIEF_OBJ),  INTENT(IN)    :: ZREF
      DOUBLE PRECISION, INTENT(INOUT) :: AC,AVAIL(NPOIN,10,NSICLA)
      DOUBLE PRECISION, INTENT(INOUT) :: ES(NPOIN,10)
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: CSTAEQ
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: FLUER_SABLE,FLUER_VASE
      DOUBLE PRECISION,  INTENT(INOUT) :: MS_SABLE(NPOIN,10) 
      DOUBLE PRECISION,  INTENT(INOUT) :: MS_VASE(NPOIN,10)
      DOUBLE PRECISION,  INTENT(INOUT) ::TOCE_MIXTE(NPOIN,10)
      DOUBLE PRECISION, INTENT(IN)      :: DT
      TYPE(BIEF_OBJ),   INTENT(IN)       ::  QSC
      INTEGER,          INTENT (IN)      :: ICQ
      DOUBLE PRECISION, INTENT(IN)     :: TOCE_VASE(10)
      !----------------------------------------------------------------!
      END SUBROUTINE 
      !================================================================!

      !================================================================!
      SUBROUTINE SUSPENSION_FREDSOE
      !----------------------------------------------------------------!
     &  (ACLADM, TAUP, NPOIN, GRAV, XMVE, XMVS, ZERO, AC,  CSTAEQ)
      USE BIEF
      TYPE(BIEF_OBJ),   INTENT(IN)    :: ACLADM, TAUP
      INTEGER,          INTENT(IN)    :: NPOIN
      DOUBLE PRECISION, INTENT(IN)    ::  GRAV, XMVE, XMVS
      DOUBLE PRECISION, INTENT(IN)    :: ZERO
      DOUBLE PRECISION, INTENT(IN)    :: AC
      TYPE(BIEF_OBJ),   INTENT(INOUT) ::  CSTAEQ
      !----------------------------------------------------------------!
      END SUBROUTINE SUSPENSION_FREDSOE
      !================================================================!


      !================================================================!
      SUBROUTINE SUSPENSION_LISTING
      !----------------------------------------------------------------!
     &  (MESH,CST,ZFCL_S,UCONV,VCONV,MASKEL,IELMT,DT,MSK,T1)
      !----------------------------------------------------------------!
      USE BIEF
      TYPE(BIEF_MESH),  INTENT(INOUT) :: MESH
      TYPE(BIEF_OBJ),   INTENT(IN)    :: CST, ZFCL_S
      TYPE(BIEF_OBJ),   INTENT(IN)    :: UCONV, VCONV, MASKEL
      INTEGER,          INTENT(IN)    :: IELMT
      DOUBLE PRECISION, INTENT(IN)    :: DT
      LOGICAL,          INTENT(IN)    :: MSK
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: T1
      !----------------------------------------------------------------!
      END SUBROUTINE SUSPENSION_LISTING
      !================================================================!


      !================================================================!
        SUBROUTINE SUSPENSION_MAIN
      ! ************************** !

     &(SLVTRA,HN,HN_TEL,MU,TOB,ACLADM,KSP,KSR,
     * CF,VOLU2D,V2DPAR,UNSV2D,AFBOR,
     & BFBOR,ZF,LICBOR,IFAMAS,MASKEL,MASKPT,U2D,V2D,NSICLA,NPOIN,
     & NPTFR,IELMT,OPTDIF,RESOL,LT,NIT,OPTBAN,OPTSUP,OPDTRA,
     & KENT,KSORT,KLOG,KINC,KNEU,KDIR,KDDL,ISOUS,NSOUS,DEBUG,
     & DTS,CSF_VASE,CSF_SABLE,ZERO,GRAV,XKX,XKY,KARMAN,
     & XMVE,XMVS,HMIN,XWC,VITCD,VITCE,PARTHENIADES,ENTET,BILMA,MSK,
     & CHARR,IMP_INFLOW_C,MESH,ZF_S,CS,CST,CTILD,CBOR,DISP,
     & IT1,IT2,IT3,IT4,TB,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,W1,
     & TE1,CLT,TE2,TE3,S,AM1_S,AM2_S,MBOR,ELAY,LIMDIF,
     & MASKTR, TETA_SUSP, AC, MASED0, MASINI, MASTEN,
     & MASTOU, ES,AVAIL,  ENTETS, PASS, ZFCL_S,
     & HPROP, FLUDPT, FLUDP, FLUER, DISP_C, KX, KY,
     & KZ, UCONV, VCONV,QSXS, QSYS, QSCLXS, QSCLYS, QSCL_S, 
     & QS_S,QS_C,CSTAEQ,ICQ,MASTCP,MASFIN,MASDEPT,MASDEP,MASSOU,
     & CORR_CONV,ZREF,SEDCO,VISC_TEL,CODE,
     & DIFT,DM1,UCONV_TEL,VCONV_TEL,ZCONV,SOLSYS,FLBOR_TEL,FLBOR_SIS,
     & FLBORTRA,NUMLIQ,NFRLIQ,MIXTE,NCOUCH_TASS,CONC_VASE,
     & TOCE_VASE,FLUER_VASE,TOCE_MIXTE,MS_SABLE,MS_VASE,TASS)
      USE BIEF_DEF
      IMPLICIT NONE
      TYPE (SLVCFG),    INTENT(INOUT) :: SLVTRA
      TYPE (BIEF_OBJ),  INTENT(IN)    :: HN,HN_TEL,MU,TOB, ACLADM
      TYPE (BIEF_OBJ),  INTENT(IN)    :: CF,VOLU2D,AFBOR,BFBOR,ZF
      TYPE (BIEF_OBJ),  INTENT(IN)    :: V2DPAR,UNSV2D,KSP,KSR
      TYPE (BIEF_OBJ),  INTENT(IN)    :: LICBOR, IFAMAS, MASKEL, MASKPT
      TYPE (BIEF_OBJ),  INTENT(IN)    :: U2D, V2D,DM1,ZCONV,FLBOR_TEL
      INTEGER,          INTENT(IN)    :: NSICLA,NPOIN,NPTFR,IELMT
      INTEGER,          INTENT(IN)    :: OPTDIF,RESOL,LT,NIT
      INTEGER,          INTENT(IN)    :: OPTBAN,OPTSUP,OPDTRA,NFRLIQ
      INTEGER,          INTENT(IN)    :: KENT, KSORT, KLOG, KINC, KNEU
      INTEGER,          INTENT(IN)    :: KDIR,KDDL,ISOUS,NSOUS
      INTEGER,          INTENT(IN)    :: DEBUG,SOLSYS,NCOUCH_TASS
      INTEGER,          INTENT(IN)    :: NUMLIQ(NFRLIQ)
      DOUBLE PRECISION, INTENT(IN)    :: DTS,CSF_VASE,CSF_SABLE,KARMAN
      DOUBLE PRECISION, INTENT(IN)    :: ZERO,GRAV,XKX,XKY
      DOUBLE PRECISION, INTENT(IN)    :: XMVE,XMVS,HMIN,XWC(NSICLA)
      DOUBLE PRECISION, INTENT(IN)    :: VITCD,VITCE,PARTHENIADES
      LOGICAL,          INTENT(IN)    :: ENTET,BILMA,MSK,CHARR
      LOGICAL,          INTENT(IN)    :: IMP_INFLOW_C
      LOGICAL,          INTENT(IN)    :: SEDCO(NSICLA),MIXTE,TASS
      TYPE(BIEF_MESH),  INTENT(INOUT) :: MESH
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: ZF_S,CS,CST,CTILD,CBOR
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: DISP,IT1,IT2,IT3,IT4,TB
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: T1,T2,T3,T4,T5,T6,T7,T8
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: T9,T10,T11,T12,W1,TE1,CLT
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: TE2,TE3,S,AM1_S,AM2_S,MBOR
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: ELAY, LIMDIF,FLBORTRA,MASKTR
      DOUBLE PRECISION, INTENT(INOUT) :: TETA_SUSP, AC(NSICLA)
      DOUBLE PRECISION, INTENT(INOUT) :: MASED0(NSICLA), MASINI(NSICLA)
      DOUBLE PRECISION, INTENT(INOUT) :: MASTEN(NSICLA), MASTOU(NSICLA)
      DOUBLE PRECISION, INTENT(INOUT) :: ES(NPOIN,10)
      DOUBLE PRECISION, INTENT(INOUT) :: TOCE_VASE(10),CONC_VASE(10)
      DOUBLE PRECISION, INTENT(INOUT) :: AVAIL(NPOIN,10,NSICLA)
      LOGICAL,          INTENT(INOUT) :: ENTETS, PASS
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: ZFCL_S,HPROP,ZREF
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: FLUDPT,FLUDP,FLUER
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: DISP_C,KX,KY,KZ,UCONV
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: VCONV,FLBOR_SIS
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: QSXS,QSYS,QSCLXS,QSCLYS
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: QSCL_S,QS_S,CSTAEQ
      TYPE (BIEF_OBJ),  INTENT(IN)    :: QS_C,VISC_TEL
      TYPE (BIEF_OBJ),  INTENT(IN)    :: UCONV_TEL,VCONV_TEL
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: FLUER_VASE,TOCE_MIXTE
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: MS_SABLE,MS_VASE
      DOUBLE PRECISION, INTENT(OUT)   :: MASTCP(NSICLA),MASFIN(NSICLA)
      DOUBLE PRECISION, INTENT(OUT)   :: MASDEPT(NSICLA),MASDEP(NSICLA)
      DOUBLE PRECISION, INTENT(OUT)   :: MASSOU
      INTEGER, INTENT(IN)             :: ICQ     
      LOGICAL, INTENT (IN)            :: CORR_CONV,DIFT
      CHARACTER(LEN=24), INTENT(IN)   :: CODE
      !----------------------------------------------------------------!
      END SUBROUTINE SUSPENSION_MAIN
      !================================================================!


      !================================================================!
      SUBROUTINE SUSPENSION_NERBED
      !----------------------------------------------------------------!
     &  (ZR, ELAY, AVA, HN, NPOIN, DT, CSF, HMIN, FLUDP, FLUER, ZF_S)
      !----------------------------------------------------------------!
      USE BIEF
      TYPE (BIEF_OBJ),  INTENT(IN)    :: ZR,ELAY,AVA,HN,ZF_S
      INTEGER,          INTENT(IN)    :: NPOIN
      DOUBLE PRECISION, INTENT(IN)    :: DT,CSF,HMIN
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: FLUDP,FLUER
      !----------------------------------------------------------------!
      END SUBROUTINE SUSPENSION_NERBED
      !================================================================!


      !================================================================!
        SUBROUTINE SUSPENSION_ROUSE 
      !================================================================!

     & (USTAR,HN,NPOIN,KARMAN,HMIN,ZERO,XWC,ZREF,T2)
      USE BIEF
      IMPLICIT NONE
      TYPE(BIEF_OBJ),   INTENT(IN) :: USTAR,HN,ZREF
      INTEGER,          INTENT(IN) :: NPOIN
      DOUBLE PRECISION, INTENT(IN) :: KARMAN,XWC,HMIN,ZERO
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: T2
      !----------------------------------------------------------------!
      END SUBROUTINE SUSPENSION_ROUSE
      !================================================================!

      !================================================================!
        SUBROUTINE SUSPENSION_BIJKER 
      !================================================================!

     &  (TAUP, HN, NPOIN, CHARR, QSC, ZREF, ZERO, HMIN, CSTAEQ, XMVE)


      USE BIEF
      IMPLICIT NONE
      TYPE(BIEF_OBJ),   INTENT(IN)    :: TAUP, HN,QSC
      TYPE (BIEF_OBJ),  INTENT(IN)    :: ZREF
      INTEGER,          INTENT(IN)    :: NPOIN
      LOGICAL,          INTENT(IN)    :: CHARR
      DOUBLE PRECISION, INTENT(IN)    :: ZERO
      DOUBLE PRECISION, INTENT(IN)    :: HMIN, XMVE
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: CSTAEQ
      !----------------------------------------------------------------!
      END SUBROUTINE SUSPENSION_BIJKER
      !================================================================!
      ! ************************** !
        SUBROUTINE SUSPENSION_EVOL
      ! ************************** !

     *  (ZFCL_S,FLUDP,FLUER,DT, NPOIN,CSF,XMVS, QFLUX,MS,
     *   SEDCO,CONC_VASE,NCOUCH_TASS)
       USE BIEF
       TYPE (BIEF_OBJ),  INTENT(INOUT) :: ZFCL_S,FLUDP,FLUER,QFLUX
       DOUBLE PRECISION, INTENT(IN)    :: DT, XMVS, CSF
       INTEGER, INTENT(IN) :: NPOIN,NCOUCH_TASS
       LOGICAL, INTENT(IN) :: SEDCO 
       DOUBLE PRECISION, INTENT(IN) :: CONC_VASE(NCOUCH_TASS)
       DOUBLE PRECISION,  INTENT(INOUT) :: MS(NPOIN,NCOUCH_TASS)
 
      END SUBROUTINE SUSPENSION_EVOL



!======================================================================!
!======================================================================!

      END INTERFACE
      END MODULE INTERFACE_SISYPHE_SUSPENSION
