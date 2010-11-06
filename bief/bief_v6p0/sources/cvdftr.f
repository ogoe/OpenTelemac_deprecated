C                       *****************
                        SUBROUTINE CVDFTR
C                       *****************
C
     *(F,FTILD,FN,FSCEXP,DIFT,ICONVF,CONV,
     * H,HN,HPROP,TETAH,UCONV,VCONV,DM1,ZCONV,SOLSYS,
     * VISC,VISC_S,SM,SMH,YASMH,SMI,YASMI,AM1,AM2,
     * ZF,FBOR,AFBOR,BFBOR,LIMTRA,MASKTR,MESH,W,TB,
     * T1,T2,T3,T4,T5,T6,T7,T10,TE1,TE2,TE3,KDIR,KDDL,KENT,DT,ENTET,
     * TETAT,AGGLOT,INFOGT,BILAN,OPTSUP,
     * ISOUSI,LT,NIT,OPDTRA,OPTBAN,MSK,MASKEL,MASKPT,MBOR,
     * S,MASSOU,OPTSOU,SLVTRA,FLBOR,V2DPAR,UNSV2D,OPTVF,FLBORTRA)
C
C
C  NOTE JMH : W IS NOT USED
C
C***********************************************************************
C  BIEF VERSION 6.0     29/12/2009    J-M HERVOUET (LNHE) 01 30 87 80 18
C                                           C MOULIN   (LNH) 30 87 83 81
C
C 27/02/2009 JMH : APPEL DE CVTFVF_POS, OPTION 14
C
C***********************************************************************
C
C  FONCTION : DIFFUSION, ADVECTION AND SOURCE TERMS FOR A TRACER
C
C  THE EQUATION SOLVED IS :
C
C
C          N+1                                            TILD
C         F           1                                  F   + DT*SM
C      ---------  -  ---  DIV ( H VISC * GRAD ( F   )) = ____________
C         DT          H                                      DT
C
C                                                      N+1
C                                  + SOURCES  + SMI * F
C                                                     ___
C                                                     H
C
C     WITH :    N+1  TILD   N
C              F   ,F     ,F  =    FONCTION DIFFUSEE
C              VISC           =    VISCOSITE TURBULENTE
C              SM             =    SECOND MEMBRE (TERMES SOURCES)
C              TETAT          =    COEFFICIENT D'IMPLICITATION
C              DT             =    PAS DE TEMPS
C                                         N+1              N
C              F              =    TETAT F  + (1-TETAT) * F
C              SMI            =    TERME SOURCE IMPLICITE.
C
C
C                   TILD      N
C     ON DISTINGUE F      ET F   AU CAS OU II Y AURAIT EU AUPARAVANT
C
C     UNE ETAPE DE PAS FRACTIONNAIRES (CONVECTION PAR EXEMPLE) DONNANT
C
C      TILD              N
C     F     A PARTIR DE F
C
C-----------------------------------------------------------------------
C
C      CONDITIONS AU LIMITES :
C
C      ==>   CONDITION DE NEUMANN
C
C      VISC DF/DN = AFBOR . F  +  BFBOR
C
C
C      ==>   CONDITION DE DIRICHLET
C
C            TRAITE PAR MODIFICATION DES EQUATIONS DANS LE SOUS
C            PROGRAMME DIRICH.
C
C-----------------------------------------------------------------------
C
C      AVERTISSEMENT :
C
C      ==>   MATDIF NE DONNE PAS LA MATRICE DE DIFFUSION, EN EFFET IL
C            MANQUE LES TERMES DE BORD ET IL Y A UN SIGNE MOINS DONT ON
C            TIENT COMPTE ICI.
C
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   F            |<-- |  F AT TIME T(N+1)
C |   FTILD        | -->|  F AFTER ADVECTION
C |   FN           | -->|  F AT TIME T(N)
C |   FSCEXP       | -->|  EXPLICIT PART OF THE SOURCE TERM
C |                |    |  EQUAL TO ZERO EVERYWHERE BUT ON SOURCES
C |                |    |  WHERE THERE IS FSCE - (1-TETAT) FN
C |                |    |  SEE DIFSOU
C |   DIFT         | -->|  IF YES, DIFFUSION IS DONE
C |   ICONVF       | -->|  OPTION FOR ADVECTION TERMS
C |                |    |  ICONVF = 1 : CHARACTERISTICS.
C |                |    |  ICONVF = 2 : S.U.P.G.
C |                |    |  ICONVF = 3 : CONSERVATIVE FINITE VOLUMES
C |                |    |  ICONVF = 4 : IDEM
C |                |    |  ICONVF = 6 : NON CONSERVATIVE PSI SCHEME.
C |                |    |  ICONVF = 7 : NON CONSERVATIVE N SCHEME.
C |                |    |  ICONVF =13 : EDGE BY EDGE FORM OF 3
C |                |    |  ICONVF =14 : IDEM
C |   CONV         | -->|  IF YES ADVECTION OF F
C |   H , HN       | -->|  DEPTH AT TIME T(N+1) AND T(N)
C |   HPROP        | -->|  WORK ARRAY
C |   TETAH        | -->|  IMPLICITATION BETWEEN H AND HN
C |   UCONV,VCONV  | -->|  ADVECTION VELOCITY FIELD
C |   DM1,ZCONV    | -->|  SEE BELOW
C |   SOLSYS       | -->|  1 OR 2. IF 2 ADVECTION FIELD IS UCONV + DM1*GRAD(ZCONV)
C |   VISC         | -->|  COEFFICIENTS DE VISCOSITE SUIVANT X,Y ET Z .
C |                |    |  SI P0 : VISCOSITE DONNEE PAR ELEMENT
C |                |    |  SINON : VISCOSITE DONNEE PAR POINT
C |   VISC_S       |<-->|  WORK ARRAY FOR SAVING VISC
C |   SM           | -->|  TERMES SOURCES .
C |   SMH          | -->|  TERME SOURCE DE L'EQUATION DE CONTINUITE
C |   YASMH        | -->|  IF YES SMH TAKEN INTO ACCOUNT
C |   SMI          | -->|  IMPLICIT SOURCE TERM
C |   YASMI        | -->|  IF YES SMI TAKEN INTO ACCOUNT
C |   AM1          |<-->|  MATRIX.
C |   AM2          |<-->|  MATRIX.
C |   ZF           | -->|  BOTTOM ELEVATION.
C |   FBOR         | -->|  CONDITIONS DE DIRICHLET SUR F.
C |   AFBOR,BFBOR  | -->|  COEFFICIENTS DES CONDITIONS DE NEUMANN
C |                |    |  VISC*DF/DN = AFBOR*F + BFBOR
C |                |    |  ATTENTION | EN P0 DONNE PAR ELEMENT DE BORD
C |                |    |              SINON DONNE PAR POINT DE BORD
C |   LIMTRA       | -->|  TYPES DE CONDITIONS AUX LIMITES SUR LES
C |                |    |  POINTS DE BORD.
C |   MASKTR(1,1)  | -->|  MASQUE VALANT 1. POUR LES SEGMENTS DIRICHLET
C |   MASKTR(1,2)  | -->|  MASQUE VALANT 1. POUR LES SEGMENTS DDL
C |   MASKTR(1,3)  | -->|  MASQUE VALANT 1. POUR LES SEGMENTS NEUMANN
C |                |    |  (ET ZERO SINON)
C |   MESH         | -->|  BLOC DES ENTIERS DU MAILLAGE.
C |   W            | -->|  TABLEAU DE TRAVAIL DE DIMENSION :
C |                |    |  NELMAX * (NOMBRE DE POINTS DANS L'ELEMENT)
C |   TB           | -->|  BLOC DE TABLEAUX DE TRAVAIL (CONTIENT T1,...)
C |   T1......T10  |<-->|  9 TABLEAUX DE TRAVAIL (PAS DE T9)
C |   TE1,TE2,TE3  |<-->|  TABLEAUX DE TRAVAIL SUR LES ELEMENTS
C |   KDIR         | -->|  CONVENTION POUR LES POINTS DE DIRICHLET
C |   KDDL         | -->|  CONVENTION POUR LES DEGRES DE LIBERTE
C |   KENT         | -->|  
C |   DT           | -->|  PAS DE TEMPS
C |   ENTET        | -->|  LOGIQUE INDIQUANT SI ON IMPRIME DES INFOS
C |                |    |  SUR LE BILAN DE MASSE DE TRACEUR
C |   TETAT        | -->|  COEFFICIENT D'IMPLICITATION DE LA CONVECTION
C |   AGGLOT       | -->|  COEFFICIENT DE MASS-LUMPING DE T.
C |   INFOGT       | -->|  LOGIQUE INDIQUANT SI ON IMPRIME DES INFOS
C |                |    |  SUR LE SOLVEUR.
C |   BILAN        | -->|  LOGIQUE INDIQUANT SI ON DOIT FAIRE UN BILAN
C |                |    |  DE MASSE. DANS CE CAS IL FAUT RETOURNER LA
C |                |    |  VALEUR DE L'APPORT DES TERMES SOURCES.
C |   OPTSUP       | <--|  VARIANTES DE LA METHODE SUPG
C |                |    |  1 : SUPG CLASSIQUE
C |                |    |  2 : SUPG MODIFIE
C |   ISOUSI       | -->|  NUMERO DE LA SOUS-ITERATION
C |   LT,NIT       | -->|  NUMERO DU PAS DE TEMPS,NOMBRE TOTAL DE PAS.
C |   OPDTRA       | -->|  MOT-CLE : OPTION POUR LA DIFFUSION DU TRACEUR
C |   OPTBAN       | -->|  OPTION DE TRAITEMENT DES BANCS DECOUVRANTS
C |                |    |  1:CLASSIQUE   2:AVEC MASQUAGE.
C |   MSK          | -->|  SI OUI, PRESENCE D'ELEMENTS MASQUES.
C |   MASKEL       | -->|  TABLEAU DE MASQUAGE DES ELEMENTS
C |                |    |  =1. : NORMAL   =0. : ELEMENT MASQUE
C |   MASKPT       | -->|  TABLEAU DE MASQUES PAR POINTS.
C |   MBOR         | -->|  MATRICE DE BORD
C |   S            | -->|  STRUCTURE BIDON
C |   MASSOU       | -->|  MASSE DE TRACEUR AJOUTEE PAR TERME SOURCE
C |                |    |  VOIR DIFSOU
C |   OPTSOU       | -->|  OPTION DE TRAITEMENT DES TERMES SOURCES.
C |                |    |  1 : NORMAL
C |                |    |  2 : DIRAC
C |   SLVTRA       | -->|  SLVCFG STRUCTURE CONTAINING DATA FOR CALLING SOLVE
C |   FLBOR        | -->|  FLUXES AT BOUNDARIES
C |   V2DPAR       | -->|  INTEGRAL OF TEST FUNCTIONS (ASSEMBLED IN PARALLEL)
C |   UNSV2D       | -->|  =1/V2DPAR
C |   FLBORTRA     |<-->|  TRACER FLUXES AT BOUNDARIES
C |   OPTVF        | -->|  OPTIONS POUR LES VOLUMES FINIS (VOIR CVTRVF)
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : 
C
C***********************************************************************
C
C   ATTENTION : POUR LES ELEMENTS DE BORD OU IL N'Y A PAS DE FROTTEMENT
C               AFBOR ET BFBOR DOIVENT ETRE NULS.
C
C   ATTENTION : A LA DISCRETISATION DE VISC
C
C**********************************************************************
C
      USE BIEF, EX_CVDFTR => CVDFTR
      USE DECLARATIONS_TELEMAC, ONLY : ADV_CAR,ADV_SUP,ADV_NSC,ADV_PSI,
     *   ADV_PSI_NC,ADV_NSC_NC,ADV_LPO,ADV_NSC_TF,ADV_PSI_TF,ADV_LPO_TF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)           :: ICONVF,ISOUSI,OPTSUP,OPDTRA,KENT
      INTEGER, INTENT(IN)           :: LT,NIT,OPTBAN,OPTSOU,KDIR,SOLSYS
      INTEGER, INTENT(IN)           :: KDDL,OPTVF
      DOUBLE PRECISION, INTENT(IN)  :: TETAT,AGGLOT,TETAH,DT
      DOUBLE PRECISION, INTENT(INOUT)  :: MASSOU
      LOGICAL, INTENT(IN)           :: INFOGT,BILAN,CONV,YASMH
      LOGICAL, INTENT(IN)           :: DIFT,MSK,ENTET,YASMI
      TYPE(SLVCFG), INTENT(INOUT)   :: SLVTRA
      TYPE(BIEF_OBJ), INTENT(IN)    :: MASKEL,MASKPT,H,HN,AFBOR,BFBOR
      TYPE(BIEF_OBJ), INTENT(INOUT) :: HPROP
      TYPE(BIEF_OBJ), INTENT(INOUT) :: F,SM,FLBORTRA
      TYPE(BIEF_OBJ), INTENT(IN)    :: FBOR,UCONV,VCONV,ZF
      TYPE(BIEF_OBJ), INTENT(IN)    :: FTILD,FN,SMI
      TYPE(BIEF_OBJ), INTENT(INOUT) :: SMH
      TYPE(BIEF_OBJ), INTENT(INOUT) :: TE1,TE2,TE3,W
      TYPE(BIEF_OBJ), INTENT(INOUT) :: T1,T2,T3,T4,T5,T6,T7,T10
      TYPE(BIEF_OBJ), INTENT(IN)    :: FSCEXP,DM1,ZCONV
      TYPE(BIEF_OBJ), INTENT(IN)    :: S,LIMTRA,FLBOR,V2DPAR,UNSV2D
      TYPE(BIEF_OBJ), INTENT(INOUT) :: VISC_S,VISC
      TYPE(BIEF_OBJ), INTENT(INOUT) :: AM1,AM2,MBOR
      TYPE(BIEF_OBJ), INTENT(INOUT) :: TB
      TYPE(BIEF_OBJ), INTENT(IN)    :: MASKTR
      TYPE(BIEF_MESH) :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION C,CFLMAX
C
      INTEGER IELMF,IELMH,IELMS,IELMU,MSKNEU,I,N,IOPT
C
      LOGICAL MSQ,FV_SCHEME    
C
      CHARACTER*16 FORMUL
C
C-----------------------------------------------------------------------
C
      IELMF = F%ELM
      IELMH = H%ELM
      IELMS = SM%ELM
      IELMU = UCONV%ELM
C
C-----------------------------------------------------------------------
C
C     DO WE HAVE A FINITE VOLUME SCHEME FOR ADVECTION ?
C
      FV_SCHEME=.FALSE.
      IF(  ICONVF.EQ.ADV_LPO.OR.ICONVF.EQ.ADV_LPO_TF.OR.
     *     ICONVF.EQ.ADV_NSC.OR.ICONVF.EQ.ADV_NSC_TF.OR.
     *     ICONVF.EQ.ADV_PSI.OR.ICONVF.EQ.ADV_PSI_TF     ) THEN
        FV_SCHEME=.TRUE.
      ENDIF
C
C-----------------------------------------------------------------------
C
C     QUAND H ET T N'ONT PAS LA MEME DISCRETISATION
C
      IF(IELMF.NE.IELMS) THEN
        CALL CHGDIS(SM ,IELMS,IELMF,MESH)
      ENDIF
      IF(IELMF.NE.IELMH.AND.YASMH) THEN
        CALL CHGDIS(SMH ,IELMH,IELMF,MESH)
      ENDIF
C
C
C-----------------------------------------------------------------------
C
C     SEMI-IMPLICITATION DE LA HAUTEUR
C     AVEC SCHEMA 5 SUR H, TETAH=0.
C
      CALL OS( 'X=CY    ' , X=HPROP , Y=H     , C=      TETAH )
      CALL OS( 'X=X+CY  ' , X=HPROP , Y=HN    , C= 1.D0-TETAH )
      CALL OS( 'X=Y     ' , X=T10   , Y=HPROP )
      IF(IELMF.NE.IELMH) THEN
        CALL CHGDIS(T10,IELMH,IELMF,MESH)
      ENDIF
C
C-----------------------------------------------------------------------
C
C     INITIALISATION DES VARIABLES
C
C     SOLUTION INITIALISEE A F AU TEMPS N
      IF(ISOUSI.EQ.1) CALL OS( 'X=Y     ' , X=F , Y=FN )
C
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C  SI SUPG, ON CONSTRUIT LA MATRICE SEMI-IMPLICITE + SUPG DANS AM2
C
      IF(ICONVF.EQ.ADV_SUP.AND.CONV) THEN
C
C       TERME EN U.GRAD(T) CENTRE :
C
        CALL MATRIX(AM2,'M=N     ','MATVGR          ',IELMF,IELMF,
     *              1.D0,S,S,S,UCONV,VCONV,S,
     *              MESH,MSK,MASKEL)
C
C       CONTRIBUTION DE SUPG AJOUTEE A AM2
C
        IF(OPTSUP.EQ.1) THEN
C         SUPG CLASSIQUE
          CALL KSUPG(TE1,TE2,1.D0,UCONV,VCONV,MESH)
          CALL MATRIX(AM2,'M=M+N   ','MASUPG          ',IELMF,IELMF,
     *                1.D0,TE1,TE2,S,UCONV,VCONV,S,MESH,MSK,MASKEL)
C
        ELSEIF(OPTSUP.EQ.2) THEN
C         SUPG MODIFIE
          CALL MATRIX(AM2,'M=M+N   ','MAUGUG          ',IELMF,IELMF,
     *                0.5D0*DT,S,S,S,UCONV,VCONV,S,MESH,MSK,MASKEL)
        ENDIF
C
C     SCHEMA N EQUATION NON CONSERVATIVE, IMPLICITE
      ELSEIF(ICONVF.EQ.ADV_NSC_NC) THEN
C
C       TERME EN U.GRAD(T) PAR SCHEMA N IMPLICITE.
C
        CALL MATRIX(AM2,'M=N     ','MATVGR         N',IELMF,IELMF,
     *              1.D0,S,S,S,UCONV,VCONV,S,MESH,MSK,MASKEL)
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C   CALCUL DE AM1 : MATRICE DE MASSE MULTIPLIEE PAR 1/DT
C
      IF(DIFT.OR..NOT.FV_SCHEME.OR..NOT.CONV.OR.BILAN) THEN
C
      IF(OPTBAN.EQ.1.OR.OPTBAN.EQ.3) THEN
        CALL OS('X=Y+Z   ',T2,ZF,HN,C)
        CALL DECVRT(TE3,T2,ZF,MESH)
C       MATRICE DE MASSE LUMPEE LOCALEMENT SUR LES BANCS DECOUVRANTS
        FORMUL='MSLUMP          '
C       MSQ SERVIRA POUR LE MASQUAGE DE LA DIFFUSION
        MSQ=.TRUE.
        IF(MSK) CALL OS('X=XY    ',TE3,MASKEL,MASKEL,C)
      ELSE
C       MATRICE DE MASSE NORMALE
        FORMUL='MATMAS          '
C       MASQUE POUR LA DIFFUSION = MASKEL
        IF(MSK) CALL OS('X=Y     ',TE3,MASKEL,MASKEL,C)
        MSQ=MSK
      ENDIF
      CALL MATRIX(AM1,'M=N     ',FORMUL,IELMF,IELMF,
     *            1.D0/DT,TE3,S,S,S,S,S,MESH,MSK,MASKEL)
C
C   MASS-LUMPING EVENTUEL
C
      IF(AGGLOT.GT.0.001D0) THEN
        CALL LUMP(T1,AM1,MESH,AGGLOT)
        CALL OM( 'M=CN    ' , AM1 , AM1 , S  , 1.D0-AGGLOT , MESH )
        CALL OM( 'M=M+D   ' , AM1 , AM1 , T1 , C           , MESH )
      ENDIF
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      IF(BILAN) THEN
C
        CALL MATVEC( 'X=AY    ',T2,AM1,SM,C,MESH)
        IF(NCSIZE.GT.1) THEN
          CALL PARCOM(T2,2,MESH)
          MASSOU = MASSOU + P_DOTS(T2,T10,MESH)
        ELSE
          MASSOU = MASSOU + DOTS(T2,T10)
        ENDIF
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C   CALCUL DES SECONDS MEMBRES
C
C     CALCUL DE DT * SM MIS DANS T2
C
C     CVTRVF AND CVTRVF_POS WILL TREAT SM IN A DIFFERENT WAY 
      IF(.NOT.FV_SCHEME) THEN
        CALL OS( 'X=CY    ' , X=T2 , Y=SM , C=DT )
      ENDIF
C
C=======================================================================
C TRAITEMENT DES DIFFERENTS TYPES DE CONVECTION :
C-----------------------------------------------------------------------
C
      IF(ICONVF.EQ.ADV_CAR.OR..NOT.CONV) THEN
C
        CALL OS( 'X=X+Y   ' , X=T2 , Y=FTILD )
        CALL MATVEC( 'X=AY    ',SM,AM1,T2,C,MESH)
C
C-----------------------------------------------------------------------
C
      ELSEIF(ICONVF.EQ.ADV_SUP.AND.CONV) THEN
C
C       AM1 RENDUE NON SYMETRIQUE SI ELLE NE L'ETAIT PAS
C
        IF(AM1%TYPEXT.NE.'Q') THEN
          CALL OM( 'M=X(M)  ' , AM1 , AM1 , S , C , MESH )
        ENDIF
C
C       CONTRIBUTION DE SUPG POUR LA MATRICE DE MASSE
C
        IF(OPTSUP.EQ.1) THEN
C         SUPG CLASSIQUE
C         TE1 ET TE2 DEJA CALCULES
          CALL MATRIX(AM1,'M=M+TN    ','MATVGR          ',IELMF,IELMF,
     *                1.D0/DT,S,S,S,TE1,TE2,S,MESH,MSK,MASKEL)
C
        ELSEIF(OPTSUP.EQ.2) THEN
C         SUPG MODIFIE
          CALL MATRIX(AM1,'M=M+TN    ','MATVGR          ',IELMF,IELMF,
     *                0.5D0,S,S,S,UCONV,VCONV,S,MESH,MSK,MASKEL)
        ENDIF
C
C       FIN DE LA CONTRIBUTION DE SUPG POUR LA MATRICE DE MASSE
C
        CALL OS( 'X=X+Y   ' , T2 , FN , FN , C )
        CALL MATVEC( 'X=AY    ',SM,AM1,T2,C,MESH)
C
C TERME DE CONVECTION EXPLICITE :
C
        CALL MATVEC( 'X=X+CAY ',SM,AM2,FN,TETAT-1.D0,MESH)
C
C ON AJOUTE A AM1 LA PARTIE DE CONVECTION IMPLICITE COMPRISE DANS AM2
C
        CALL OM( 'M=M+CN  ' , AM1,AM2 , S , TETAT , MESH )
C
C-----------------------------------------------------------------------
C
      ELSEIF(ICONVF.EQ.ADV_NSC_NC.AND.CONV) THEN
C
C       AM1 RENDUE NON SYMETRIQUE SI ELLE NE L'ETAIT PAS
C
        IF(AM1%TYPEXT.NE.'Q') THEN
          CALL OM( 'M=X(M)  ' , AM1 , AM1 , S , C , MESH )
        ENDIF
C
        CALL OS( 'X=X+Y   ' , T2 , FN , FN , C )
        CALL MATVEC( 'X=AY    ',SM,AM1,T2,C,MESH)
C
C TERME DE CONVECTION EXPLICITE :
C
        CALL MATVEC( 'X=X+CAY ',SM,AM2,FN,TETAT-1.D0,MESH)
C
C ON AJOUTE A AM1 LA PARTIE DE CONVECTION IMPLICITE COMPRISE DANS AM2
C
        CALL OM( 'M=M+CN  ' , AM1,AM2 , S , TETAT , MESH )
C
C-----------------------------------------------------------------------
C
      ELSEIF(ICONVF.EQ.ADV_PSI_NC.AND.CONV) THEN
C
C SCHEMA PSI
C
C       TERME AM1 * FN CLASSIQUE
C
        CALL OS( 'X=X+Y   ' , T2 , FN , FN , C )
        CALL MATVEC( 'X=AY    ',SM,AM1,T2,C,MESH)
C
C       TERME DE CONVECTION EXPLICITE PAR SCHEMA PSI
C
        CALL VGFPSI(T5,IELMF,UCONV,VCONV,FN,DT,-1.D0,CFLMAX,
     *              T6,T7,MESH,MSK,MASKEL)
        CALL OS( 'X=X+Y   ' , SM , T5 , T5 , C )
C
C-----------------------------------------------------------------------
C
      ELSEIF( (ICONVF.EQ.ADV_LPO.OR.
     *         ICONVF.EQ.ADV_NSC.OR.
     *         ICONVF.EQ.ADV_PSI    ).AND.CONV ) THEN
C
C CONSERVATIVE EQUATION, DISTRIBUTIVE SCHEMES (LEO POSTMA, N AND PSI)
C                        LEO POSTMA AND N-SCHEME ARE THE SAME IN 2D
C
C       TO BE REMOVED WHEN ALL CALLS TO CVDFTR ARE CHECKED
C       OPTVF SHOULD BE 0 (VELOCITY FIELD OBEYS THE CONTINUITY EQUATION)
C       OR 10 (VELOCITY FIELD DOES NOT OBEY THE CONTINUITY EQUATION)
        IOPT=10*(OPTVF/10)
C       OPTION FOR DISTRIBUTING THE FLUXES (HERE 2 OR 3)
        IF(ICONVF.EQ.ADV_LPO) IOPT=IOPT+2
        IF(ICONVF.EQ.ADV_NSC) IOPT=IOPT+2
        IF(ICONVF.EQ.ADV_PSI) IOPT=IOPT+3
C
        IF(TB%N.LT.22) THEN
          WRITE(LU,*) 'SIZE OF TB TOO SMALL IN CVDFTR'
          CALL PLANTE(1)
          STOP
        ENDIF
        CALL CVTRVF(F,FN,FSCEXP,DIFT,CONV,H,HN,HPROP,UCONV,VCONV,
     *              DM1,ZCONV,SOLSYS,VISC,VISC_S,SM,SMH,YASMH,SMI,YASMI,
     *              FBOR,MASKTR,MESH,
     *              TB%ADR(13)%P,TB%ADR(14)%P,TB%ADR(15)%P,
     *              TB%ADR(16)%P,TB%ADR(17)%P,TB%ADR(18)%P,
     *              TB%ADR(19)%P,TB%ADR(20)%P,TB%ADR(21)%P,
     *              TB%ADR(22)%P,
     *              AGGLOT,TE1,DT,ENTET,BILAN,
     *              OPDTRA,MSK,MASKEL,S,MASSOU,OPTSOU,
C                                                       YAFLBOR
     *              LIMTRA%I,KDIR,KDDL,MESH%NPTFR,FLBOR,.TRUE.,
     *              V2DPAR,UNSV2D,IOPT,FLBORTRA,MASKPT)
C       SI ON REPART ICI ON NE FAIT PAS LES DIRICHLET, CA MARCHE AUSSI
C       ET ON PEUT VERIFIER EXACTEMENT LA CONSERVATION DE LA MASSE
        IF(.NOT.DIFT) RETURN
        CALL MATVEC( 'X=AY    ',SM,AM1,F,C,MESH)
C
C-----------------------------------------------------------------------
C
      ELSEIF( (ICONVF.EQ.ADV_LPO_TF.OR.
     *         ICONVF.EQ.ADV_NSC_TF.OR.
     *         ICONVF.EQ.ADV_PSI_TF    ).AND.CONV ) THEN
C
C EDGE-BASED VERSIONS, FOR TIDAL FLATS
C CONSERVATIVE EQUATION, DISTRIBUTIVE SCHEMES (LEO POSTMA, N AND PSI)
C                        LEO POSTMA AND N-SCHEME ARE THE SAME IN 2D
C
C       OPTION FOR DISTRIBUTING THE FLUXES (HERE 2 OR 3)
        IF(ICONVF.EQ.ADV_LPO_TF) IOPT=2
        IF(ICONVF.EQ.ADV_NSC_TF) IOPT=2
        IF(ICONVF.EQ.ADV_PSI_TF) IOPT=3
        IF(TB%N.LT.22) THEN
          WRITE(LU,*) 'SIZE OF TB TOO SMALL IN CVDFTR'
          CALL PLANTE(1)
          STOP
        ENDIF
        CALL CVTRVF_POS(F,FN,FSCEXP,DIFT,CONV,H,HN,HPROP,UCONV,VCONV,
     *              DM1,ZCONV,SOLSYS,VISC,VISC_S,SM,SMH,YASMH,SMI,YASMI,
     *              FBOR,MASKTR,MESH,
     *              TB%ADR(13)%P,TB%ADR(14)%P,TB%ADR(15)%P,
     *              TB%ADR(16)%P,TB%ADR(17)%P,TB%ADR(18)%P,
     *              TB%ADR(19)%P,TB%ADR(20)%P,TB%ADR(21)%P,
     *              TB%ADR(22)%P,
     *              AGGLOT,TE1,DT,ENTET,BILAN,
     *              OPDTRA,MSK,MASKEL,S,MASSOU,OPTSOU,
C                                                       YAFLBOR
     *              LIMTRA%I,KDIR,KDDL,MESH%NPTFR,FLBOR,.TRUE.,
     *              V2DPAR,UNSV2D,IOPT,FLBORTRA,MASKPT,
     *            MESH%GLOSEG%I(                 1:  MESH%GLOSEG%DIM1),
     *            MESH%GLOSEG%I(MESH%GLOSEG%DIM1+1:2*MESH%GLOSEG%DIM1),
     *            MESH%NBOR%I)
C       SI ON REPART ICI ON NE FAIT PAS LES DIRICHLET, CA MARCHE AUSSI
C       ET ON PEUT VERIFIER EXACTEMENT LA CONSERVATION DE LA MASSE
        IF(.NOT.DIFT) RETURN
        CALL MATVEC( 'X=AY    ',SM,AM1,F,C,MESH)
      ELSE
C
C-----------------------------------------------------------------------
C
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'CVDFTR : OPTION DE CONVECTION INCONNUE : ',ICONVF
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) 'CVDFTR: UNKNOWN ADVECTION OPTION : ',ICONVF
        ENDIF
        CALL PLANTE(1)
        STOP      
C
C-----------------------------------------------------------------------
C
      ENDIF
C
C FIN DU TRAITEMENT DES DIFFERENTS TYPES DE CONVECTION :
C=======================================================================
C
C                   CALCUL DES MATRICES
C
C-----------------------------------------------------------------------
C   CALCUL DE AM2 : - MATRICE DE DIFFUSION , ET TERMES DE BORD
C
      IF(DIFT) THEN
C
        IF(OPDTRA.EQ.2) THEN
C             SAUVEGARDE DE LA DIFFUSION
              CALL OS('X=Y     ',VISC_S,VISC,VISC,C)
C             MULTIPLICATION PAR HPROP DE LA DIFFUSION
           CALL OV_2('X=XY    ',VISC%R,1,T10%R,1,T10%R,1,C,
     *                          VISC%MAXDIM1,VISC%DIM1)
           IF(VISC%DIM2.EQ.3) THEN
           CALL OV_2('X=XY    ',VISC%R,2,T10%R,1,T10%R,1,C,
     *                          VISC%MAXDIM1,VISC%DIM1)
           CALL OV_2('X=XY    ',VISC%R,3,T10%R,1,T10%R,1,C,
     *                          VISC%MAXDIM1,VISC%DIM1)
           ENDIF
        ENDIF
C
C       CALCUL DE LA MATRICE DE DIFFUSION (OPTION AVEC MONOTONIE)
C
        CALL MATRIX(AM2,'M=N     ','MATDIF       MON',IELMF,IELMF,
     *              1.D0,S,S,S,VISC,S,S,MESH,MSQ,TE3)
C
        IF(OPDTRA.EQ.2) THEN
C         MULTIPLICATION DE LA MATRICE PAR 1/HPROP
          CALL OS( 'X=1/Y   ',T4,T10,T10,C,
     *             IOPT=2,INFINI=0.D0,ZERO=1.D-2)
          CALL OM( 'M=X(M)  ' , AM2 , AM2 , S  , C , MESH )
          CALL OM( 'M=DM    ' , AM2 , AM2 , T4 , C , MESH )
C         RESTITUTION DE LA DIFFUSION
          CALL OS('X=Y     ',VISC,VISC_S,VISC_S,C)
        ENDIF
C
C   PRISE EN COMPTE DES TERMES DE BORD DANS LA MATRICE DE DIFFUSION
C
        MSKNEU=3
        CALL MATRIX(MBOR,'M=N     ','FMATMA          ',
     *              IELBOR(IELMF,1),IELBOR(IELMF,1),
     *              -1.D0,AFBOR,S,S,S,S,S,
     *              MESH,.TRUE.,MASKTR%ADR(MSKNEU)%P)
        CALL OM( 'M=M+N   ' , AM2 , MBOR , S , C , MESH )
C
C       TERME DE DIFFUSION EXPLICITE
C
        CALL MATVEC( 'X=AY    ',T1,AM2,FN,C,MESH)
        CALL OS( 'X=X+CY  ' , SM , T1 , T1 , TETAT-1.D0 )
C
C       TERME DE DIFFUSION IMPLICITE ( AM1 + TETAT * AM2 )
C
        IF(AM1%TYPEXT.NE.'Q'.AND.AM2%TYPEXT.EQ.'Q') THEN
          CALL OM( 'M=X(M)  ' , AM1 , AM1 , S , C , MESH )
        ENDIF
        CALL OM( 'M=M+CN  ' , AM1,AM2 , S , TETAT , MESH )
C
C       TERMES DE CONTRAINTE AU BORD
C
        CALL VECTOR(T2,'=','MASVEC          ',IELBOR(IELMF,1),
     *              1.D0,BFBOR,S,S,S,S,S,MESH,
     *              .TRUE.,MASKTR%ADR(MSKNEU)%P)
        CALL OSDB( 'X=X+Y   ' , SM , T2 , T2 , C , MESH )
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C   PRISE EN COMPTE DU TERME IMPLICITE VENANT DES SOURCES PONCTUELLES :
C
      IF(YASMH.AND..NOT.(FV_SCHEME.AND.CONV)) THEN
C
        IF(OPTSOU.EQ.1) THEN
C         MODIF JMH DU 23/09/98
          CALL VECTOR(T2,'=','MASVEC          ',IELMF,
     *                1.D0,SMH,S,S,S,S,S,MESH,MSK,MASKEL)
          CALL OS( 'X=Y/Z   ' ,T1,T2,T10,C,
     *              IOPT=2,INFINI=0.D0,ZERO=1.D-3)
C         PARTIE IMPLICITE DU TERME SOURCE PONCTUEL                       
C         - TETAT T 1/HPROP SOMME ( SCE PSI D(OMEGA)                      
          CALL OS( 'X=CX    ' , T1 , T1 , T1 , TETAT )                  
          CALL OM( 'M=M+D   ' , AM1 , AM1 , T1 , TETAT , MESH )
C         PREPARATION DE LA PARTIE EXPLICITE 
          CALL OS( 'X=YZ    ' , T1 , SMH , FSCEXP , C )
          CALL VECTOR(T2,'=','MASVEC          ',IELMF,
     *                1.D0,T1,S,S,S,S,S,MESH,MSK,MASKEL)
          CALL OS( 'X=Y/Z   ' ,T1,T2,T10,C,
     *             IOPT=2,INFINI=0.D0,ZERO=1.D-3)
          CALL OS( 'X=X+Y   ' , SM , T1 , T1 , C )             
        ELSEIF(OPTSOU.EQ.2) THEN
          CALL OS( 'X=Y/Z   ' ,T1,SMH,T10,C,
     *              IOPT=2,INFINI=0.D0,ZERO=1.D-3)
C         PARTIE EXPLICITE DU TERME SOURCE PONCTUEL
C         1/HPROP (FSCE-(1-TETAT)FN) SMH
          CALL OS( 'X=X+YZ  ' , SM , T1 , FSCEXP , C )
C         PARTIE IMPLICITE DU TERME SOURCE PONCTUEL
C         - TETAT T 1/HPROP SOMME ( SCE PSI D(OMEGA)
          CALL OS( 'X=CX    ' , T1 , T1 , T1 , TETAT )
          CALL OM( 'M=M+D   ' , AM1 , AM1 , T1 , TETAT , MESH )
        ENDIF
C
      ENDIF
C
C   IMPLICITE TERM IF THERE IS ONE :
C
C   THE TREATMENT BELOW ENSURES THAT IF THE EXPLICIT SOURCE TERM
C   IS IN THE FORM  K*FN/H AND SMI EQUALS -K THEN THE TWO TERMS
C   WILL BE BALANCED (CASE OF EROSION AND DEPOSITION).
C
C                  FV_SCHEME : IMPLICIT SOURCE TERM HAS BEEN TREATED 
C                  AND WILL NOT BE DONE TWICE
      IF(YASMI.AND..NOT.(FV_SCHEME.AND.CONV)) THEN
        CALL MATRIX(AM2,'M=N     ','MATMAS          ',IELMF,IELMF,
     *              -1.D0,S,S,S,S,S,S,MESH,MSK,MASKEL)
C       MASS-LUMPING EVENTUEL
        IF(AGGLOT.GT.0.001D0) THEN
          CALL LUMP(T1,AM2,MESH,AGGLOT)
          CALL OM( 'M=CN    ' , AM2 , AM2 , S  , 1.D0-AGGLOT , MESH )
          CALL OM( 'M=M+D   ' , AM2 , AM2 , T1 , C           , MESH )
        ENDIF
C       CALCUL DE SMI/H (TAILLE SUFFISANTE DE H NON VERIFIEE!!)
        IF(OPTBAN.GT.0) THEN
C         DIVISION BY H WITH HARDCODED CLIPPING AT 0.01
          CALL CPSTVC(SMI,T4)
          DO I=1,SMI%DIM1
            IF(T10%R(I).LT.1.D-2) THEN
              T4%R(I)=0.D0
            ELSE
              T4%R(I)=SMI%R(I)/H%R(I)
            ENDIF
          ENDDO
        ELSE
C         DIVISION WITHOUT CHECKING
          CALL OS( 'X=Y/Z   ',X=T4,Y=SMI,Z=H)
        ENDIF
        CALL OM( 'M=X(M)  ' , AM2 , AM2 , S  , C , MESH )
        CALL OM( 'M=MD    ' , AM2 , AM2 , T4 , C , MESH )
C       AJOUT A LA MATRICE AM1
        IF(AM1%TYPEXT.NE.'Q') THEN
          CALL OM( 'M=X(M)  ' , AM1 , AM1 , S , C , MESH )
        ENDIF
        CALL OM( 'M=M+N   ' , AM1 , AM2 , S , C , MESH )
      ENDIF
C
C-----------------------------------------------------------------------
C
      IF(ICONVF.EQ.ADV_CAR.AND..NOT.DIFT) THEN
        CALL OS( 'X=Y     ' , F , FTILD , FTILD , C )
      ENDIF
      IF(ICONVF.EQ.ADV_PSI_NC.AND.CONV) THEN
        CALL LUMP(T1,AM1,MESH,1.D0)
        IF(NCSIZE.GT.1) THEN
          CALL PARCOM(T1,2,MESH)
          CALL OS( 'X=Y     ' , T2 , SM , SM , C )
          CALL PARCOM(T2,2,MESH)
          CALL OS( 'X=Y/Z   ' , F , T2 , T1 , C ,2,0.D0,1.D-6)
        ELSE
          CALL OS( 'X=Y/Z   ' , F , SM , T1 , C ,2,0.D0,1.D-6)
        ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
C
C   CONDITIONS AUX LIMITES AUX BORDS (POINTS DE TYPE DIRICHLET).
C
      CALL DIRICH(F, AM1, SM,FBOR,LIMTRA%I,TB,MESH,KDIR,MSK,MASKPT)
C
C-----------------------------------------------------------------------
C
C   RESOLUTION DU SYSTEME LINEAIRE :
C
      CALL SOLVE(F,AM1,SM,TB,SLVTRA,INFOGT,MESH,AM2)
C
C-----------------------------------------------------------------------
C
C     CALCUL DU FLUX DE TRACEUR AU BORD
C
      IF(.NOT.FV_SCHEME) THEN
        DO I=1,MESH%NPTFR
          N=MESH%NBOR%I(I)
          FLBORTRA%R(I)=FLBOR%R(I)*(TETAT*F%R(N)+(1.D0-TETAT)*FN%R(N))
        ENDDO
!     ELSE
!       FLBORTRA ALREADY COMPUTED    
      ENDIF
C
C-----------------------------------------------------------------------
C
      IF(IELMF.NE.IELMS) CALL CHGDIS(SM  ,IELMF,IELMS,MESH)
      IF(IELMF.NE.IELMH.AND.YASMH) CALL CHGDIS(SMH ,IELMF,IELMH,MESH)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
