!! modifs 05/04/2010
! cstaeq tient compte du % de la couche qq soit la formule utilisee
      ! ********************************* !
        SUBROUTINE SUSPENSION_COMPUTATION
      ! ********************************* !
!
     &  (SLVTRA, HN,HN_TEL,UCONV, VCONV,  CF, MU,TOB,ACLADM, KSP,KSR,
     &   ELAY, AVA, AFBOR, BFBOR, LIMDIF, CLT, MASKEL, MASKTR,
     &   MASKPT, IFAMAS, NPOIN, IELMT, NPTFR, ITRA, LT, NIT, RESOL,
     &   OPTBAN, KENT,KDDL,KDIR,KSORT,KLOG,KINC,KNEU,
     &   OPTSUP, OPDTRA, DEBUG, CSF_VASE,CSF_SABLE,
     &   TETA_SUSP, DT, MASED0, ZERO, XWC, KARMAN, XMVE, XMVS, GRAV,
     &   HMIN, VITCD, VITCE,PARTHENIADES, ENTETS,
     &   ENTET,BILMA,MSK,CHARR,IMP_INFLOW_C,MESH,ZF,CS,
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

C**********************************************************************
C SISYPHE VERSION 6.0  22/12/04  F. HUVELIN                            
C
C 05/05/2008 : ADAPTATION AU CONVECTEUR VOLUMES FINIS
C 09/05/2008 : FLUDP ENLEVE DE SUSPENSION_FLUX, 
C              SUSPENSION_NERBED SUPPRIME.
C 28/05/2008 : NOUVEAU SUSPENSION_BILAN AVEC FLUX PAR FRONTIERES
C 09/06/2008 : NOUVEAU SUSPENSION_BILAN AVEC FLBORTRA FOURNI PAR CVDFTR
C 12/06/2008 : INVERSION DES SECTIONS "TREATING SMALL DEPTHS" ET
C                                     "LIMITATION OF FLUER..."
C 25/06/2008 : APPEL A DIFFIN QUI ETAIT DANS SUSPENSION_MAIN
C 31/07/2008 : APPEL A SUSPENSION_FLUX CASSE EN DEUX : DEPOT + EROSION
C 16/09/2009 : AVAIL(NPOIN,10,NSICLA)
C**********************************************************************
!
!
!               ! ====================================== !
!               ! Main subroutine for the computation of !
!               ! the concentration and the elevation    !
!               ! Solving the equation :                 !
!               !    D(ZF)                               !
!               !    ----  + DIV(QS) + (E-D)Za = 0       !
!               !    DT                                  !
!               ! if coupling, DIV(QS) already computed  !
!               ! else         DIV(QS) = 0               !
!               ! ====================================== !
!
!
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
C |   HPROP        | <= | Work array for the water depth               C
C |   DISP_C       | <= | Work array for the dispersion                C
C |________________|____|______________________________________________C
C                    <=  Can't be change by the user                   C
C                    =>  Can be changed by the user                    C 
C ---------------------------------------------------------------------C
!                                                                      !
! CALLED BY SUSPENSION_MAIN                                            !
!                                                                      !
! CALL      SUSPENSION_ZREF
!           SUSPENSION_CONV
!           CHARAC (BIEF)                                              !
!           DIFFCL (BIEF)                                              !
!           SUSPENSION_FLUX                                            !
!           SUSPENSION_NERBED                                          !
!           CVDFTR (BIEF)                                              !
!           SUSPENSION_LISTING                                         !
!           SUSPENSION_BILAN                                           !
!                                                                      !
!======================================================================!
!======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS               !
!======================================================================!
!======================================================================!

      ! 1/ MODULES
      ! ----------
      USE INTERFACE_SISYPHE, 
     &    EX_SUSPENSION_COMPUTATION => SUSPENSION_COMPUTATION
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU

      ! 2/ GLOBAL VARIABLES
      ! -------------------
      TYPE (SLVCFG),    INTENT(INOUT) :: SLVTRA
      TYPE (BIEF_OBJ),  INTENT(IN)    :: ZF,VOLU2D,V2DPAR,UNSV2D
      TYPE (BIEF_OBJ),  INTENT(IN), TARGET    :: HN,HN_TEL
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: UCONV,VCONV
      TYPE (BIEF_OBJ),  INTENT(IN)    :: MU,KSP,KSR
      TYPE (BIEF_OBJ),  INTENT(IN)    :: CF,TOB,ACLADM,ELAY,LICBOR
      TYPE (BIEF_OBJ),  INTENT(IN)    :: AFBOR,BFBOR
      TYPE (BIEF_OBJ),  INTENT(IN)    :: MASKEL,MASKPT,IFAMAS
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: MASKTR,LIMDIF,CLT
      INTEGER,          INTENT(IN)    :: NPOIN,IELMT,NPTFR,ITRA,LT
      INTEGER,          INTENT(IN)    :: NIT,RESOL,OPTBAN,KENT,KDDL
      INTEGER,          INTENT(IN)    :: KDIR,OPTSUP,OPDTRA,SOLSYS
      INTEGER,          INTENT(IN)    :: KSORT,KLOG,KINC,KNEU
      INTEGER,          INTENT(IN)    :: NFRLIQ,NSICLA,NCOUCH_TASS
      INTEGER,          INTENT(IN)    :: DEBUG
      INTEGER,          INTENT(IN)    :: NUMLIQ(NFRLIQ)
      DOUBLE PRECISION, INTENT(IN)    :: CSF_VASE, TETA_SUSP, DT, MASED0
      DOUBLE PRECISION, INTENT(IN)    :: ZERO, XWC, CSF_SABLE,AVA(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: KARMAN, XMVE, XMVS, GRAV, HMIN
      DOUBLE PRECISION, INTENT(IN)    :: VITCD,VITCE,PARTHENIADES
      LOGICAL,          INTENT(IN)    :: ENTETS,ENTET,BILMA,MSK,SEDCO
      LOGICAL,          INTENT(IN)    :: CHARR, IMP_INFLOW_C,CORR_CONV
      LOGICAL,          INTENT(IN)    :: DIFT,MIXTE, TASS
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
C
!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!
!
! 3/ LOCAL VARIABLES
! ------------------
! 
      INTEGER          :: I,K,SOLSYS_SIS,OPTVF,BID(1),RESOL_MOD
      DOUBLE PRECISION :: TETAH,AGGLOT
      DOUBLE PRECISION :: CSF
      LOGICAL YASMI2
      TYPE (BIEF_OBJ),  POINTER :: HOLD
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: SAVE_UCONV,SAVE_VCONV
      SAVE_UCONV=>UCONV%R
      SAVE_VCONV=>VCONV%R
!
!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!
!
!     JMH 18/04/2008 :
!     PRISE EN COMPTE DES DETAILS DE L'EQUATION DE CONTINUITE
!     DANS TELEMAC-2D OU 3D AVEC SOLSYS=2 ON UTILISE DM1 ET ZCONV
!     QUI NE SONT PAS CONSTRUITS QUAND ON APPELLE SISYPHE
!     AU PREMIER PAS DE TEMPS
!     DE PLUS UCONV_TEL ET VCONV_TEL SONT ALORS DIFFERENTS
!     DE U2D ET V2D
!
      SOLSYS_SIS=1
      IF(LT.NE.1.AND.CODE(1:7).EQ.'TELEMAC') SOLSYS_SIS=SOLSYS
!     VOIR PLUS BAS UNE MODIFICATION DE SOLSYS_SIS : IF(CORR_CONV)...
!
!
!
! TWO OPTIONS : ICQ=1: FREDSOE REFERENCE CONC. ZREF = 2.D50
!               ICQ=2: BIJKER METHOD ZREF= MAX (KSP,KS)        
!       
! 
!     COMPUTATION OF THE REFERENCE ELEVATION  -->  ZREF
!
        IF(ICQ.EQ.1) CALL OS('X=Y     ', X=ZREF, Y=KSP)
        IF(ICQ.EQ.2) CALL OS('X=Y     ', X=ZREF, Y=KSR)

!
!     COMPUTATION OF THE ADVECTION VELOCITIES
!     PRISE EN COMPTE DU PROFIL VERTICAL DES CONCENTRATIONS ET DES VITESSES
! 
!     OPTVF : CHIFFRE DES DIZAINES  0 : NORMAL
!                                   1 : CHAMP CONVECTEUR NE VERIFIE PAS LA
!                                       CONTINUITE
!
!     OPTVF : CHIFFRE DES UNITES    0 : CONSTANTE NULLE
!                                   1 : CONSTANTE CHI-TUAN
!                                   2 : CONSTANTE LEO POSTMA
!                                   VOIR CVTRVF DANS BIEF ET
!                                   RELEASE NOTES 5.7
!
!
      IF(CORR_CONV.AND.(.NOT.SEDCO)) THEN
! 
!        ICI ON PREND PROVISOIREMENT U2D ET V2D 
!        MEME SI SOLSYS_TEL=2 CAR IL FAUDRAIT CORRIGER AUSSI DM1
         SOLSYS_SIS=1
         CALL CPSTVC(U2D,T12)        
         CALL SUSPENSION_CONV(TOB,XMVE, KSR,NPOIN,ZREF,U2D,V2D,HN,HMIN,
     *                        UCONV,VCONV,KARMAN,ZERO,XWC,T1,T12)
!        FORME DU CONVECTEUR QUI ACCEPTE UN CHAMP CONVECTEUR
!        QUI N'OBEIT PAS A LA CONTINUITE + CONSTANTE CHI-TUAN
         OPTVF=11
!
      ELSE
!
!       ICI ON UTILISE DES POINTEURS POUR EVITER LA COPIE
!
        IF(SOLSYS_SIS.EQ.1) THEN
          UCONV%R=>U2D%R
          VCONV%R=>V2D%R
        ELSE
!         DANS CE CAS UCONV_TEL A ETE TRANSMIS
          UCONV%R=>UCONV_TEL%R
          VCONV%R=>VCONV_TEL%R
        ENDIF
!        FORME DU CONVECTEUR QUI VEUT UN CHAMP CONVECTEUR
!        QUI OBEIT A LA CONTINUITE + CONSTANTE CHI-TUAN
         OPTVF=1
!
      ENDIF
!
!     ADVECTION WITH CHARACTERISTICS
!      
      IF(RESOL == 1) THEN
         IF (DEBUG > 0) WRITE(LU,*) 'CHARAC'
         CALL CHARAC(CS,CTILD,1,UCONV,VCONV,S,S,DT,IFAMAS,
     &          IELMT, NPOIN,1,1,MSK,MASKEL,AM1_S%X,AM1_S%D,
     &               TB,IT1%I,IT2%I,IT3%I,IT4%I,
     &               MESH,MESH%NELEM,MESH%NELMAX,MESH%IKLE,
     &               MESH%SURDET)
         IF (DEBUG > 0) WRITE(LU,*) 'END_CHARAC'
      ENDIF
!
!     FLUX COMPUTATION AT THE TIME N (CALLED FOR THE FIRST CLASS
!                                     WHICH MUST BE THE SAND) 
! 
! le dépôt est traité de la même manière pour les mélanges et les 
!  sédiments simples
!!---> T2: rapppport entre coententration au fond et concentration moy
!

      IF (DEBUG > 0) WRITE(LU,*) 'SUSPENSION_DEPOT'
      CALL SUSPENSION_DEPOT(TOB, HN,ACLADM,NPOIN, HMIN,XWC,
     &  VITCD, ZERO,KARMAN,XMVE,T1,T2,ZREF,FLUDPT,DEBUG,SEDCO)     
!     
! +++++++++++++++++++++++++++++
! L'érosion est traitée différemment   : passer tass en argument  
!  FROTT DE PEAU TAUP  --> T4
! 
      CALL OS('X=CYZ   ', X= T4, Y= TOB, Z= MU, C=1.D0)
      CALL OS('X=+(Y,C)', X=T4, Y=T4, C=ZERO)
! V6P0 changements CV
      IF(.NOT.MIXTE) THEN
        IF(.NOT.SEDCO) THEN       
         IF (DEBUG > 0) WRITE(LU,*) 'SUSPENSION_EROSION'
          CALL SUSPENSION_EROSION(T4,HN,ACLADM,AVA, 
     &     NPOIN,CHARR,XMVE,XMVS,GRAV,HMIN,XWC,ZERO,
     &     ZREF,AC,FLUER,CSTAEQ,QS_C,ICQ,DEBUG)
         IF (DEBUG > 0) WRITE(LU,*) 'END_SUSPENSION_EROSION'
!
! this should be included in SUSPENSION_EROSION
! 
         DO I=1,NPOIN
           FLUER%R(I)=MIN(FLUER%R(I),ELAY%R(I)*AVA(I)/DT*CSF_SABLE)  
         ENDDO
!
! vase pure : new subroutine
! 
        ELSE
          CALL  SUSPENSION_EROSION_COH (T4,NPOIN,
     &       XMVE,XMVS,GRAV, VITCE, PARTHENIADES,ZERO, DEBUG, 
     &      FLUER, ES, TOCE_VASE, NCOUCH_TASS, DT, MS_VASE%R,TASS)
!
          IF(.NOT.TASS) THEN          
            DO I=1,NPOIN
              FLUER%R(I)=MIN(FLUER%R(I),ELAY%R(I)*CSF_VASE/DT)
            ENDDO
          ENDIF 
        ENDIF
! SEDIMENT mixte 
       ELSE
!         CSF=CSF_VASE ! defini dans sediment mixte
         IF(.NOT.SEDCO) THEN
          IF(DEBUG > 0) WRITE(LU,*) 'SUSPENSION_FLUX_MIXTE'
          CALL SUSPENSION_FLUX_MIXTE(T4,HN,ACLADM,CS, 
     &         NPOIN,CHARR,XMVE,XMVS,GRAV,HMIN,XWC,
     &         ZERO,KARMAN,PARTHENIADES,FLUER,
     &         FLUER_VASE,ZREF,AC,CSTAEQ,QS_C,ICQ,DEBUG,
     &         AVAIL,NSICLA,ES,TOCE_VASE,NCOUCH_TASS,DT,
     &         TOCE_MIXTE%R,MS_SABLE%R,MS_VASE%R)
          IF (DEBUG > 0) WRITE(LU,*) 'END_SUSPENSION_FLUX_MOY'
        ENDIF
        IF(SEDCO) CALL OS('X=Y     ',X=FLUER, Y=FLUER_VASE)         
      ENDIF
!
! ....fin modif sediments mixtes
!
!     TREATING SMALL DEPTHS      
!
      IF(OPTBAN.EQ.1) THEN
        DO I = 1, NPOIN
          IF(HN%R(I).LE.HMIN) THEN
            FLUDPT%R(I)=0.D0
            FLUER%R(I) =0.D0
          ENDIF       
        ENDDO
      ELSEIF(OPTBAN.EQ.2) THEN
        CALL OS('X=XY    ',X=FLUER ,Y=MASKPT)
        CALL OS('X=XY    ',X=FLUDPT,Y=MASKPT)
      ENDIF
!
!     IMPLICIT SOURCE TERM FOR THE DEPOSITION 
!
      CALL OS('X=-Y    ',X=T9,Y=FLUDPT)
!
!     EXPLICIT SOURCE TERM WITHOUT PUNCTUAL SOURCES 
!
      CALL OS('X=Y     ',X=T11,Y=FLUER)
!
      DO I=1,NPOIN
        IF(HN%R(I).GT.HMIN) THEN
          T11%R(I)=T11%R(I)/HN%R(I)
        ELSE
          T11%R(I)=0.D0
        ENDIF
      ENDDO
!
!     JMH: 01/08/2005
!     IN DIFFIN A SPECIFIC TREATMENT IS DONE IF THE ADVECTION METHOD
!     IS THE CHARACTERISTICS: FREE OUTPUTS ARE TREATED LIKE DIRICHLET
!     THIS SPECIFIC TREATMENT IS CANCELLED HERE BY SENDING A MODIFIED
!     VALUE OF RESOL : RESOL_MOD (IN DIFFIN THE ONLY TEST IS:
!     IF(RESOL.EQ.1) THEN .... ELSE ....  ENDIF)
!
      RESOL_MOD=RESOL
      IF(RESOL_MOD.EQ.1) RESOL_MOD=2
      CALL DIFFIN(MASKTR,LIMDIF%I,LICBOR%I,CLT%I,U2D%R,V2D%R,
     &            MESH%XNEBOR%R,MESH%YNEBOR%R,
     &            MESH%NBOR%I,MESH%KP1BOR%I,NPTFR,
     &            KENT,KSORT,KLOG,KINC,KNEU,KDIR,KDDL,RESOL_MOD,
     &            MESH%NELBOR%I,NPOIN,MESH%NELMAX,
     &            MSK,MASKEL%R,0,
     &            .FALSE.,BID,BID,BID,CS,CBOR,MESH)
!
!     IMPOSE THE EQUIL. CONCENTRATION FOR THE INFLOW NODES  ! 
!     HERE CBOR FROM BOUNDARY CONDITIONS FILE OR SUBROUTINE CONLIT
!     OVERWRITTEN
!
!     T2 = rapport entre conc au fond et conc moy. doit etre conservé jusqu'à ce stade
!
      IF(IMP_INFLOW_C) THEN
!
        DO K = 1, NPTFR
          IF(CLT%I(K).EQ.KENT) THEN
            I = MESH%NBOR%I(K)
            IF(.NOT.SEDCO) THEN
CV voir modif dans suspension_erosion
CV               CBOR%R(K) = CSTAEQ%R(I)/T2%R(I)*AVA(I)
               CBOR%R(K) = CSTAEQ%R(I)/T2%R(I)
               IF(MIXTE) CBOR%R(K) = FLUER%R(I)/T2%R(I)/XWC
            ELSE
               CBOR%R(K) = FLUER%R(I)/XWC
            ENDIF
          ENDIF
        ENDDO
!
      ENDIF
!
!
!     ADVECTION-DISPERSION STEP
!
!     CONFIGURATION FOR CALLING CVDFTR
!
      TETAH  = 1.D0 - TETA_SUSP
      MASSOU = 0.D0
      AGGLOT=1.D0
      YASMI2 = .TRUE.
!
!     LES FLUX AUX BORDS DOIVENT ETRE DONNES A CVDFTR (CAS DES VOLUMES FINIS)
!     ET A SUSPENSION_BILAN
!     AVEC SISYPHE SEUL : IL FAUT LES CALCULER
!     EN COUPLAGE       : FOURNI PAR LES CODES APPELANTS SAUF A LA PREMIERE
!                         ITERATION
!
      IF(CODE(1:7).NE.'TELEMAC'.OR.LT.EQ.1) THEN
        CALL VECTOR(FLBOR_SIS,'=','FLUBDF          ',IELBOR(IELMT,1),
!                        HPROP (ICI HPROP=HN, A VOIR)
     *              1.D0,HN   ,HN,HN,UCONV,VCONV,VCONV,
     *              MESH,.FALSE.,MASKPT)          
      ELSE
        CALL OS('X=Y     ',X=FLBOR_SIS,Y=FLBOR_TEL)
!       ON DOIT CHANGER AUSSI LES FLUX AU BORD SI ON A
!       UNE CORRECTION DU CONVECTEUR (T12 DOIT AVOIR ETE
!       CONSERVE DEPUIS L'APPEL A SUSPENSION_CONV)
        IF(CORR_CONV.AND.(.NOT.SEDCO)) THEN
          CALL OSBD('X=CXY   ',FLBOR_SIS,T12,T12,1.D0,MESH)
        ENDIF
      ENDIF       
!
!     LE VRAI H DU PAS PRECEDENT EST UTILISE PAR LE CONVECTEUR VOLUMES FINIS
      IF(CODE(1:7).EQ.'TELEMAC') THEN
        IF(OPTBAN.NE.0) THEN
          CALL CPSTVC(CST,T1)
!         HN_TEL N'EST PAS CLIPPE
          DO I=1,NPOIN
            T1%R(I)=MAX(HN_TEL%R(I),HMIN)
          ENDDO
          HOLD=>T1
        ELSE
          HOLD=>HN_TEL
        ENDIF
      ELSE
!       DANS CE CAS H ET HN SONT CONFONDUS
        HOLD=>HN
      ENDIF
!
      IF (DEBUG > 0) WRITE(LU,*) 'APPEL DE CVDFTR'
      CALL CVDFTR 
     & (CST, CTILD, CS, T2,
!                            H         HTILD
     &  DIFT, RESOL, .TRUE., HN, HOLD, HPROP, TETAH, 
     &  UCONV,VCONV,DM1,ZCONV,SOLSYS_SIS,
!                     TEXP SMH  YASMH   TIMP
     &  DISP, DISP_C, T11, T2, .FALSE., T9,  YASMI2,AM1_S,AM2_S,
     &  ZF, CBOR, AFBOR, BFBOR, LIMDIF, MASKTR, MESH,
     &  W1, TB, T8, T12, T3, T4, T5, T6, T7, T10, TE1, TE2, TE3,
!                                                       BILMAS
     &  KDIR,KDDL,KENT,DT,ENTET,TETA_SUSP,AGGLOT,ENTET,.FALSE.,OPTSUP,
     &  1, LT, NIT, OPDTRA, OPTBAN, MSK, MASKEL, MASKPT, MBOR, S,
!               OPTSOU
     &  MASSOU, 1,     SLVTRA,FLBOR_SIS,V2DPAR,UNSV2D,OPTVF,FLBORTRA)
      IF (DEBUG > 0) WRITE(LU,*) 'END_CVDFTR'
!
      DO I=1,NPOIN
        FLUDP%R(I)=FLUDPT%R(I)*CST%R(I)
      ENDDO
!
!     EVOLUTION COMPUTATION AND UPDATING DATA 
! PASSER TASS EN ARGUMENT
!
!      IF(.NOT.MIXTE.AND..NOT.TASS) THEN    
!        CALL OS('X=Y-Z   ', X=ZFCL_S, Y=FLUDP, Z=FLUER)
!        CALL OS('X=CX    ', X=ZFCL_S, C=DT/CSF)
!      ELSE
!        CALL SIS_ERODE(ZFCL_S, FLUDP, FLUER,DT,
!     *                 NPOIN,XMVS,XKV, T3,SEDCO,
!     *                 CONC_VASE,NCOUCH_TASS,
!     *                 MS_SABLE%R,MS_VASE%R)     
!      ENDIF
       IF(.NOT.SEDCO) THEN
          CALL SUSPENSION_EVOL(ZFCL_S, FLUDP, FLUER,DT,
     *                 NPOIN,CSF_SABLE, XMVS,T3,MS_SABLE%R,
     *                 SEDCO,CONC_VASE,NCOUCH_TASS) 
       ELSE 
          CALL SUSPENSION_EVOL(ZFCL_S, FLUDP, FLUER,DT,
     *                 NPOIN,CSF_VASE,XMVS, T3,MS_VASE%R,
     *                 SEDCO,CONC_VASE,NCOUCH_TASS) 
       ENDIF        
!
!     LISTING OF THE MIN/MAX VALUES 
!     
      IF (ENTETS) THEN
         IF (DEBUG > 0) WRITE(LU,*) 'SUSPENSION_LISTING'
         CALL SUSPENSION_LISTING
     &        (MESH,CST,ZFCL_S,UCONV,VCONV,MASKEL,
     &         IELMT,DT,MSK,T1)
         IF (DEBUG > 0) WRITE(LU,*) 'END_SUSPENSION_LISTING'
      ENDIF
!
!     MASS-BALANCE FOR THE SUSPENSION 
!      
      IF(BILMA) THEN
         IF (SEDCO) CSF = CSF_VASE
         IF(.NOT.SEDCO) CSF = CSF_SABLE
         IF (DEBUG > 0) WRITE(LU,*) 'SUSPENSION_BILAN'
         CALL SUSPENSION_BILAN
     &        (MESH,CST,HN,ZFCL_S,MASKEL,IELMT,ITRA,LT,NIT,
     &         DT,CSF,MASSOU,MASED0,MSK,ENTETS,MASTEN,MASTOU,
     &         MASINI,T2,T3,MASFIN,MASDEPT,MASDEP,AGGLOT,VOLU2D,
     &         NUMLIQ,NFRLIQ,NPTFR,FLBORTRA)
         IF (DEBUG > 0) WRITE(LU,*) 'END_SUSPENSION_BILAN'
      ENDIF
!
      UCONV%R=>SAVE_UCONV
      VCONV%R=>SAVE_VCONV
!      
!======================================================================!
!======================================================================!
!
      RETURN
      END
