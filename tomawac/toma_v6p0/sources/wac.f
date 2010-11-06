C                         **************
                          SUBROUTINE WAC
C                         **************
C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  EDF/DER - LABORATOIRE NATIONAL D'HYDRAULIQUE  -  CHATOU (FRANCE)   C
C                                                                     C
C  CODE TOMAWAC DE MODELISATION DES ETATS DE MER EN ELEMENTS FINIS    C
C                                                                     C
C  TOMAWAC : RESOLUTION DE L'EQUATION D'EVOLUTION DE LA DENSITE       C
C  ********  SPECTRO-ANGULAIRE DE VARIANCE (A PARTIR DE LA CONSERVA-  C
C            TION DE L'ACTION D'ONDE DANS LE CAS GENERAL)             C
C                                                                     C
C                       ---- VERSION 5.2 ----                         C
C                                                                     C
C   - CREE POUR VERSION 1.0  LE 01/02/95 PAR F. MARCOS                C
C   - MOD. POUR VERSION 1.2  LE 07/11/96 PAR M. BENOIT                C
C   - PASSAGE FORTRAN 90     LE 25/08/00 PAR OPTIMER (02 98 44 24 51) C
C     POUR VERSION 5.1       SUIVI EDF/LNHE ASSURE PAR D. VIOLEAU     C
C   - FUSION TOMAWAC/COWADIS LE 12/06/01 PAR OPTIMER (02 98 44 24 51) C
C     POUR VERSION 5.2       SUIVI EDF/LNHE ASSURE PAR M. BENOIT      C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  APPELS :    - PROGRAMME(S) APPELANT  :  HOMERE_TOMAWAC             C
C  ********    - PROGRAMME(S) APPELE(S) :  PREDA2, OV    , CORFON,    C
C                                          LECSUI, CONDIW, VARTEL,    C
C                                          PRENL1, USTAR1, USTAR2,    C
C                                          INITAB, IMPR  , LIMWAC,    C
C                                          INIPHY, PRELEO, DUMP2D,    C
C                                          DUMP3D, PREPRO, PROPA ,    C
C                                          SEMIMP, TRANSF, SOR3D      C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  REMARQUES :                                                        C
C  ***********                                                        C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C
C----------------------------------------------------------------------
C     DECLARARATIONS
C----------------------------------------------------------------------
C
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_TOMAWAC
C
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/ LNG,LU
C
C     VARIABLES DECLAREES LOCALEMENT DANS LA PROCEDURE.
C
      INTEGER          NVARCL, LT , LT1 , NRK   , NPV   , NPC , NPM,
CMB...Modif MB/GM nb points de sortie      
     *                 NOLEO(99)
      LOGICAL          ISLEO(99)
CMB...fin de modif      
      INTEGER          DATE(3)       , TIME(3) 
      INTEGER          NPTL  , IBID, IP
C
C     TV1 TEMPS CORRESPONDANT AU VENT 1
C     TV2 TEMPS CORRESPONDANT AU VENT 2
C
      DOUBLE PRECISION LAMBD0, C     , Z(1)  , DEUPI , DTSI,
     *                 AT    , TV1   , TV2   , TC1 , TC2 , TM1 , TM2
      DOUBLE PRECISION VITVEN, VITMIN
      INTEGER, ALLOCATABLE :: QINDI(:)
      INTEGER          NBD , ADC , MDC , JDC , HDC, I1, I2, NVHMA,NVCOU
      LOGICAL          IMPRES, DEBRES
C
C
C
C=====C
C  1  C INITIALISATIONS DES VARIABLES LOCALES
C=====C======================================
C
C.....1.2 CONSTANTES PHYSIQUES.
C     """""""""""""""""""""""""
      GRAVIT=9.806D0
      DEUPI =2.D0*3.14159265358979D0
C
C.....1.3 PARAMETRES GENERAUX DU CALCUL.
C     """""""""""""""""""""""""""""""""
      NVARCL=0
      NPTL  =NPOIN3*NF
C
C.....1.4 INITIALISATION DES TABLEAUX DATE ET TIME
C     """"""""""""""""""""""""""""""""""""""""""""
      ADC=INT(DDC*1.D-8)
      MDC=INT(DDC*1.D-6)
      JDC=INT(DDC*1.D-4)
      HDC=INT(DDC*1.D-2)
      DATE(1)=ADC
      DATE(2)=MDC-100*ADC
      DATE(3)=JDC-100*MDC
      TIME(1)=HDC-100*JDC
      TIME(2)=INT(DDC-100.D0*HDC)
      TIME(3)=0
C
C=====C
C  2  C LECTURE DES CONDITIONS LIMITES ET INDICES DES POINTS FRONTIERES
C=====C================================================================
C
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE LECLIM'
      CALL LECLIM_TOMAWAC 
     *(SLIFBR%I,SFBOR%R,NPTFR,NBOR,STDGEO,WAC_FILES(WACCLI)%LU,
!BD_INCKA modif //
!    *  MESH%ISEG%I,MESH%XSEG%R,MESH%YSEG%R,MESH%NACHB%I )
     *  MESH%ISEG%I,MESH%XSEG%R,MESH%YSEG%R,MESH%NACHB%I ,MESH,
     *  BOUNDARY_COLOUR%I)
      IF(DEBUG.GT.0) WRITE(LU,*) 'SORTIE DE LECLIM'
!BD_INCKA fin modif //
C
C-----------------------------------------------------------------------
C
C COMPLEMENT DE LA STRUCTURE DE DONNEES POUR BIEF
C
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE INBIEF'
      CALL INBIEF(SLIFBR%I,KLOG,SITR31,SITR32,SITR33,
     *            LVMAC,IELM2,LAMBD0,SPHE,MESH,ST1,ST2,1,1,EQUA)
      IF(DEBUG.GT.0) WRITE(LU,*) 'SORTIE DE INBIEF'
C
C-----------------------------------------------------------------------
C
C LECTURE DE LA COTE DU FOND (ZF) SUR LE FICHIER DE GEOMETRIE
C
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE FONSTR'
      CALL FONSTR (ST1,SZF,ST2,ST3,WAC_FILES(WACGEO)%LU,
     *             WAC_FILES(WACFON)%LU,WAC_FILES(WACFON)%NAME,MESH,
     *             1.D0,.TRUE.)
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE FONSTR'
C
C CORRECTION EVENTUELLE DES VALEURS DU FOND (OU CALCUL DU FOND SI CELA
C N'A PAS ETE FAIT DANS FONSTR)
C EN STANDARD, CORFON NE FAIT RIEN (ATTENTION, ALLER CHERCHER LE CORFON
C DE TOMAWAC)
C
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE CORFON'
      CALL CORFON
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE CORFON'
C
C.....CALCUL DE LA PROFONDEUR D'EAU (TABLEAU DEPTH)
C
      DO IP=1,NPOIN2
        DEPTH(IP)=ZREPOS-ZF(IP)
        IF(DEPTH(IP).LT.PROMIN) DEPTH(IP)=0.9D0*PROMIN
      ENDDO
C
C-----------------------------------------------------------------------
C
C PREPARATION DES SORTIES GRAPHIQUES
C
        ! CREATION DU JEU DE DONNEES POUR UN FORMAT DE FICHIER
        ! FORMAT_RES.
        ! LE JEU DE DONNEES EST CREE DANS LE FICHIER NRES, ET EST
        ! DEFINIT PAR UN TITRE ET LES VARIABLES A ECRIRE. 
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE CREATE_DATASET'
        CALL CREATE_DATASET(WAC_FILES(WACRES)%FMT, ! FORMAT FICHIER RESULTAT
     *                      WAC_FILES(WACRES)%LU,  ! LU FICHIER RESULTAT
     *                      TITCAS,     ! TITRE DE L'ETUDE
     *                      MAXVAR,     ! MAX VARIABLES SORTIE
     *                      TEXTE,      ! NOMS VARIABLES SORTIE
     *                      SORLEO)     ! SORTIE OU PAS DES VARIABLES
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE CREATE_DATASET'
        ! ECRITURE DU MAILLAGE DANS LE FICHIER SORTIE :
        ! SI ON EST ON PARALLEL, FAUT L'INDIQUER VIA NCSIZE ET NPTIR.
        ! LES AUTRES INFORMATIONS SONT DANS MESH.
        ! EN PLUS : DATE/TEMPS DE DEPART ET LES COORDONNEES DE
        ! L'ORIGINE.
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE WRITE_MESH'
        CALL WRITE_MESH(WAC_FILES(WACRES)%FMT, ! FORMAT FICHIER RESULTAT     
     *                  WAC_FILES(WACRES)%LU,  ! LU FICHIER RESULTAT
     *                  MESH,          ! DESCRIPTEUR MAILLAGE
     *                  1,             ! NOMBRE DE PLAN /NA/
     *                  DATE,          ! DATE DEBUT
     *                  TIME,          ! HEURE DEBUT
     *                  I_ORIG,J_ORIG) ! COORDONNEES DE L'ORIGINE.
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE WRITE_MESH'
C
C-----------------------------------------------------------------------
C
C CONDITIONS INITIALES
C
C
C
C=====C INITIALISATION DES VECTEURS DE DISCRETISATION, DU COURANT,
C  2  C DU VENT ET DU SPECTRE DE VARIANCE.
C=====C===========================================================
      LT=0
      DTSI=DT/NSITS
C
      IF (SUIT) THEN
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE LECSUI'
        CALL LECSUI
     *( SF%R    , NPLAN   , NF      , STETA%R, SFR%R  ,
     *  NELEM2  , NPOIN2  , AT      , SUC%R  , SVC%R  ,
     *  SUC1%R  , SVC1%R  , SUC2%R  , SVC2%R , SUV%R  ,
     *  SVV%R   , SUV1%R  , SVV1%R  , SUV2%R , SVV2%R ,
     *  VENT    , TV1     , TV2     , COUSTA , 
     *  WAC_FILES(WACPRE)%LU ,
     *  BINPRE  , STRA02%R, SDEPTH%R, TC1 , TC2 , ZM1 , ZM2 ,
     *  SDZHDT%R, TM1     , TM2     , MAREE  )
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE LECSUI'
      ELSE
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE CONDIW'
        CALL CONDIW
     *( AT, LT , DEUPI , TC1 , TC2 , NPC, TV1, TV2, NPV, TM1, TM2, NPM,
     *  NVHMA  , NVCOU )
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE CONDIW'
      ENDIF
C
C
C
C=====C
C  3  C UTILISATION (EVENTUELLE) DE LA VARIABLE TELEMAC.
C=====C=================================================
      IF (DONTEL) THEN
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE VARTEL'
        CALL VARTEL
     *( STRA31%R, MESH%X%R, MESH%Y%R, SDEPTH%R,
     *  SUC%R   , SVC%R   , ZREPOS     , STRA32%R,
     *  SF%R    , NPLAN      , NF         , NPOIN2     )
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE VARTEL'
      ENDIF
C
      DO IP=1,NPOIN2
         IF(DEPTH(IP).LE.0.D0) THEN
           IF(LNG.EQ.1) THEN
             WRITE(LU,*) ''
             WRITE(LU,*) '*************************'
             WRITE(LU,*) ' ! PROFONDEUR NEGATIVE ! '
             WRITE(LU,*) '   ARRET DU PROGRAMME    '
             WRITE(LU,*) '*************************'
             CALL PLANTE(0)
           ELSE
             WRITE(LU,*) ''
             WRITE(LU,*) '**************************'
             WRITE(LU,*) ' ! NEGATIVE WATER DEPTH ! '
             WRITE(LU,*) '   END OF THE COMPUTATION '
             WRITE(LU,*) '**************************'
             CALL PLANTE(0)
           ENDIF
         ENDIF
      ENDDO
C
C
C
C=====C
C  4  C CALCULS PREPARATOIRES POUR INTERACTIONS NON-LINEAIRES.
C=====C=======================================================
      IF(STRIF.EQ.1) THEN
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE PRENL1'
        CALL PRENL1( IANGNL, COEFNL, NPLAN , NF , RAISF , XLAMD )
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE PRENL1'
      ENDIF
C
      ALLOCATE(QINDI(NPLAN))
      IF (STRIA.EQ.2) THEN
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE PREQT2'
        CALL PREQT2(STETA%R,NPLAN,BDISPB,BDSSPB,NBD,QINDI)
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE PREQT2'
      ENDIF
C
C=====C INITIALISATION DE LA CONTRAINTE DE HOULE, PUIS CALCUL DES
C  5  C VITESSE DE FROTTEMENT U*, RUGOSITE Z0 ET DIRECTION INITIALES.
C=====C==============================================================
C
C.....5.1 INITIALISATION DE LA CONTRAINTE DE HOULE INITIALE.
C     """"""""""""""""""""""""""""""""""""""""""""""""""""""
      CALL OV ( 'X=C     ' ,STRA41%R, STRA32%R , STRA33%R,
     *                      0.D0 , NPOIN2 )
C
C.....5.2 CALCUL DE U* ET Z0 SELON LA METHODE CONSIDEREE.
C     """""""""""""""""""""""""""""""""""""""""""""""""""
      IF (VENT) THEN
        IF (SVENT.EQ.1) THEN
          IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE USTAR1'
          CALL USTAR1
     *( STRA42%R , STRA44%R , STRA41%R , SUV%R , SVV%R ,
     *  CDRAG    , ALPHA    , XKAPPA   , ZVENT , GRAVIT,
     *  NPOIN2   )
          IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE USTAR1'
        ELSEIF (SVENT.EQ.2) THEN
          IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE USTAR2'
          CALL USTAR2(STRA42%R,SUV%R,SVV%R,NPOIN2)
          IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE USTAR2'
        ELSE
          IF (LNG.EQ.1) THEN
            WRITE(LU,*)
     *      'PB DANS WAC : VENT PRESENT, MAIS SVENT NON CORRECT'
          ELSE
            WRITE(LU,*)
     *      'PB IN WAC : WIND PRESENT, BUT SVENT NOT CORRECT'
          ENDIF
          CALL PLANTE(1)
          STOP
        ENDIF
      ENDIF
C
C.....5.3 CALCUL DE LA DIRECTION DU VENT
C     """"""""""""""""""""""""""""""""""
      VITMIN=1.D-3
      IF (VENT) THEN
        DO IP=1,NPOIN2
          VITVEN=SQRT(SUV%R(IP)**2+SVV%R(IP)**2)
          IF (VITVEN.GT.VITMIN) THEN
            STRA43%R(IP)=ATAN2(SUV%R(IP),SVV%R(IP))
          ELSE
            STRA43%R(IP)=0.D0
          ENDIF
        ENDDO
      ENDIF
C
C
C
C=====C
C  6  C INITIALISATION DE CERTAINS TABLEAUX UTILES.
C=====C============================================
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE INITAB'
      CALL INITAB( SIBOR%I, MESH%IFABOR%I, NELEM2)
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE INITAB'
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE IMPR'
      CALL IMPR(LISPRD,LT,AT,LT,3)
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE IMPR'
C
C
C
C=====C
C  7  C AFFECTATION DES CONDITIONS AUX LIMITES A L'INSTANT INITIAL.
C=====C============================================================
      I1=NF+1
      I2=NF+NPLAN
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE LIMWAC'
      CALL LIMWAC
     *(SF%R    , SFBOR%R , SLIFBR%I , NPTFR  , NPLAN , NF    ,
     * STETA%R , SFR%R   , NPOIN2   , NBOR   , AT    , LT    ,
     * DDC     , LIMSPE  , FPMAXL   , FETCHL , SIGMAL, SIGMBL,
     * GAMMAL  , FPICL   , HM0L     , APHILL , TETA1L, SPRE1L,
     * TETA2L  , SPRE2L  , XLAMDL   , MESH%X%R   ,
     * MESH%Y%R, KENT    , KSORT    , WAC_FILES(WACFO1)%LU  , 
     * WAC_FILES(WACBI1)%LU         , BINBI1 ,
     * SUV%R   , SVV%R   , SPEULI   , VENT  , VENSTA, GRAVIT ,
     * DEUPI   , SPRIVE%R, NPRIV    , STRA31%R , STRA32%R,
     * DEPTH   , FRABL   ,BOUNDARY_COLOUR%I)
       IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE LIMWAC'
C
C
C
C=====C CALCUL DES NOMBRES D'ONDE (XK), DE LA VITESSE DE GROUPE (CG) ET
C  8  C DU FACTEUR DE PASSAGE (B) EN SPECTRE DE VARIANCE EN (FR,TETA).
C=====C=================================================================
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE INIPHY'
      CALL INIPHY
     *( SXK%R   , SCG%R , SB%R , SDEPTH%R , SFR%R ,
     *  SCOSF%R , NPOIN2   , NF      , PROINF      , SPHE     )
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE INIPHY'
C
C=====C
C  8b C MISE A ZERO DU SPECTRE SUR LES POINTS OU PROF < PROMIN
C=====C=======================================================
C
      IF(.NOT.PROINF) THEN
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE ECRETE'
        CALL ECRETE(SF%R,SDEPTH%R,NPOIN2,NPLAN,NF,PROMIN)
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE ECRETE'
      ENDIF
C
C
C=====C
C  9  C SORTIES GRAPHIQUES (EVENTUELLES) A L'ETAT INITIAL.
C=====C===================================================
C
C.....9.1 CHOIX DES POINTS DE SORTIE DU SPECTRE DIRECTIONNEL.
C     """""""""""""""""""""""""""""""""""""""""""""""""""""""
C
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE PRELEO_MPI'
      CALL PRELEO_MPI
     *(XLEO,YLEO,NPLEO,MESH%X%R,MESH%Y%R,MESH%IKLE%I,MESH%SURDET%R,
     * NPOIN2,NELEM2,NOLEO,ISLEO)
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE PRELEO_MPI'
C
C.....9.2 TEST POUR SAVOIR SI ON IMPRIME OU PAS.
C     """"""""""""""""""""""""""""""""""""""""""
      IMPRES=.FALSE.
      DEBRES=.FALSE.
      IF (LT.EQ.GRADEB) THEN
        IMPRES=.TRUE.
        DEBRES=.TRUE.
      ENDIF
C
      IF (IMPRES) THEN
C
C.....9.3 IMPRESSION (EVENTUELLE) DES VARIABLES SUR LE MAILLAGE 2D.
C     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
C
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE DUMP2D'
        CALL DUMP2D(LT,DEUPI,SF%R,NPTL)
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE DUMP2D'
C
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE BIEF_DESIMP' 
        CALL BIEF_DESIMP(WAC_FILES(WACRES)%FMT,VARSOR,
     *            HIST,0,NPOIN2,WAC_FILES(WACRES)%LU,'STD',AT,
     *            LT,GRAPRD,GRAPRD,
     *            SORLEO,SORIMP,MAXVAR,TEXTE,GRADEB,GRADEB)
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE BIEF_DESIMP' 
C
C.....9.4 IMPRESSION (EVENTUELLE) DES SPECTRES DIRECTIONNELS.
C     """""""""""""""""""""""""""""""""""""""""""""""""""""""
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE ECRSPE' 
        CALL ECRSPE
     *( SF%R    , STRA01%R , STETA%R, NPLAN ,
     *  SFR%R   , NF  , NF   , NPOIN2      , AT ,LT ,
     *  STRA01%R, SITR01%I, NOLEO , NPLEO , WAC_FILES(WACLEO)%LU ,
     *  WAC_FILES(WACLEO)%FMT, DEBRES , TITCAS , DATE , TIME ,
     *  ISLEO , MESH%KNOLG%I )
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE ECRSPE'
C
      ENDIF
C
C
C=====C
C  10 C PREPARATION DE LA PROPAGATION (REMONTEE DES CARACTERISTIQUES).
C=====C===============================================================
      IF (PROP) THEN
        CALL IMPR(LISPRD,LT,AT,LT,1)
        CALL IMPR(LISPRD,LT,AT,LT,2)
        NRK=3
C
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE PREPRO'         
        CALL PREPRO
     * ( SCX%R    , SCY%R     , SCT%R    , SCF%R   , DT    ,
     *   NRK      , MESH%X%R  , MESH%Y%R , STETA%R ,
     *   SCOSTE%R , SSINTE%R  , SFR%R    , MESH%IKLE%I     ,
     *   SIBOR%I  , SETAP1%I  , STRA01%R , SSHP1%R ,
     *   SSHP2%R  , SSHP3%R   , SSHZ%R   , SSHF%R  ,
     *   SELT%I   , SETA%I    , SFRE%I   , SDEPTH%R,
     *   SDZHDT%R , SDZX%R    , SDZY%R   , SUC%R   ,
     *   SVC%R    , SDUX%R    , SDUY%R   , SDVX%R  ,
     *   SDVY%R   , SXK%R     , SCG%R    , SCOSF%R ,
     *   STGF%R   , SITR01%I  , NPOIN3   , NPOIN2  , NELEM2,
     *   NPLAN    , NF        , MESH%SURDET%R, COURAN,
     *   SPHE     , PROINF   , PROMIN,MESH)
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE PREPRO' 
      ENDIF
C
C=====C
C  11 C BOUCLE EN TEMPS PRINCIPALE.
C=====C============================
      DO 10 LT=1,NIT
C
C.....11.1 AFFECTATION DE LA DATE DE FIN DU PAS DE TEMPS COURANT.
C     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""
      AT=AT+DT
      CALL IMPR(LISPRD,LT,AT,LT,3)
C
C.....11.2 AFFECTATION DES CONDITIONS AUX LIMITES.
C     """"""""""""""""""""""""""""""""""""""""""""
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE LIMWAC' 
      CALL LIMWAC
     *(SF%R    , SFBOR%R , SLIFBR%I , NPTFR  , NPLAN , NF    ,
     * STETA%R , SFR%R   , NPOIN2      , NBOR   , AT    , LT    ,
     * DDC        , LIMSPE     , FPMAXL      , FETCHL , SIGMAL, SIGMBL,
     * GAMMAL     , FPICL      , HM0L        , APHILL , TETA1L, SPRE1L,
     * TETA2L     , SPRE2L     , XLAMDL      , MESH%X%R   ,
     * MESH%Y%R   , KENT       , KSORT       , WAC_FILES(WACFO1)%LU , 
     * WAC_FILES(WACBI1)%LU    , BINBI1 ,
     * SUV%R   , SVV%R   , SPEULI      , VENT  , VENSTA, GRAVIT ,
     * DEUPI      , SPRIVE%R, NPRIV      , STRA31%R , STRA32%R,
     * DEPTH      , FRABL      ,BOUNDARY_COLOUR%I)
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE LIMWAC'
C
C.....11.2b MISE A ZERO DU SPECTRE SUR LES POINTS OU PROF < PROMIN
C     """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""
      IF (.NOT.PROINF) THEN
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE ECRETE'
        CALL ECRETE( SF%R    , SDEPTH%R, NPOIN2, NPLAN , NF , PROMIN)
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE ECRETE'
      ENDIF
C
      IF (MAREE) THEN
       LT1=MAX((LT/LAM)*LAM,2)
       IF (LT.EQ.LT1) THEN
        DO IP=1,NPOIN2
          DEPTH(IP)=ZREPOS-ZF(IP)
        ENDDO
C
C......11.3 A JOUR DE LA BATHY ET DES COURANTS
C      """"""""""""""""""""""""""""""""""""""""
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE CORMAR'       
        CALL CORMAR
     *( AT    , LT    , TC1   , TC2   , TV1   , TV2   , TM1   , TM2   ,
     *  NPC   , NPM   , NVHMA , NVCOU )
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE CORMAR'
        DO IP=1,NPOIN2
          IF (DEPTH(IP).LT.PROMIN) DEPTH(IP)=0.9D0*PROMIN
        ENDDO
C
C......11.4 PREPARATION DE LA PROPAGATION (REMONTEE DES CARACTERISTIQUES).
C      """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE INIPHY'
        CALL INIPHY
     * ( SXK%R  , SCG%R , SB%R , SDEPTH%R , SFR%R  ,
     *   SCOSF%R, NPOIN2   , NF      , PROINF      , SPHE      )
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE INIPHY'
C
        IF (PROP) THEN
         CALL IMPR(LISPRD,LT,AT,LT,1)
         CALL IMPR(LISPRD,LT,AT,LT,2)
         NRK=3
         IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE PREPRO'
         CALL PREPRO
     * ( SCX%R    , SCY%R     , SCT%R    , SCF%R   , DT    ,
     *   NRK      , MESH%X%R  , MESH%Y%R , STETA%R ,
     *   SCOSTE%R , SSINTE%R  , SFR%R    , MESH%IKLE%I     ,
     *   SIBOR%I  , SETAP1%I  , STRA01%R , SSHP1%R ,
     *   SSHP2%R  , SSHP3%R   , SSHZ%R   , SSHF%R  ,
     *   SELT%I   , SETA%I    , SFRE%I   , SDEPTH%R,
     *   SDZHDT%R , SDZX%R    , SDZY%R   , SUC%R   ,
     *   SVC%R    , SDUX%R    , SDUY%R   , SDVX%R  ,
     *   SDVY%R   , SXK%R     , SCG%R    , SCOSF%R ,
     *   STGF%R   , SITR01%I  , NPOIN3   , NPOIN2  , NELEM2,
     *   NPLAN    , NF        , MESH%SURDET%R, COURAN,
     *   SPHE     , PROINF    , PROMIN   , MESH)
         IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE PREPRO'
        ENDIF
       ENDIF
      ENDIF
C
C.....11.3 PROPAGATION (INTERPOLATION AU PIED DES CARACTERISTIQUES).
C     """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
      IF (PROP) THEN
        CALL IMPR(LISPRD,LT,AT,LT,5)
! 
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE PROPA'  
        CALL PROPA
     *( SF%R       , SB%R    , SSHP1%R, SSHP2%R, SSHP3%R,
     *  SSHZ%R     , SSHF%R  , SELT%I , SETA%I , SFRE%I ,
     *  MESH%IKLE%I, SETAP1%I, NPOIN3    , NPOIN2    , NELEM2,
     *  NPLAN         , NF   , COURAN    ,STRA01%R   , STRA02%R )
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE PROPA'
       ENDIF

C
C.....11.4 INTEGRATION DES TERMES SOURCES.
C     """"""""""""""""""""""""""""""""""""
      IF (TSOU) THEN
        CALL IMPR(LISPRD,LT,AT,NSITS,4)
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE SEMIMP'
        CALL SEMIMP(SF%R,SXK%R,SFR%R,SDFR%R,SDEPTH%R,SUV%R,SVV%R ,
     *  MESH%X%R      , MESH%Y%R      , 
     *  WAC_FILES(WACVEB)%LU  , WAC_FILES(WACVEF)%LU  , NBOR  , NPTFR ,
     *  DDC   , TV1   , TV2   , NPV   , SXRELV%R      , SYRELV%R      ,
     *  SUV1%R, SVV1%R, SUV2%R, SVV2%R, STETA%R       , SSINTE%R      ,
     *  SCOSTE%R      , INDIV , TAILF , RAISF , GRAVIT, CFROT1, CMOUT1,
     *  CMOUT2, AT    , DTSI  , ROAIR , ROEAU , XKAPPA, BETAM , DECAL ,
     *  CDRAG , ALPHA , ZVENT , NF    , NPLAN , NPOIN2, IANGNL, COEFNL,
     *  F1    , NSITS , SMOUT , SFROT , SVENT , STRIF , VENT  , VENSTA,
     *  VX_CTE, VY_CTE, SBREK , ALFABJ, GAMBJ1, GAMBJ2, IQBBJ , IHMBJ ,
     *  IFRBJ , BORETG, GAMATG, IWHTG , IFRTG , ALFARO, GAMARO, GAM2RO,
     *  IDISRO, IEXPRO, IFRRO , BETAIH, EM2SIH, IFRIH , COEFHS, XDTBRK,
     *  NDTBRK, STRIA , ALFLTA, RFMLTA, KSPB  , BDISPB, BDSSPB, PROINF,
     *  DF_LIM, LIMIT , CIMPLI,
     *  WAC_FILES(WACVEB)%NAME, WAC_FILES(WACVEF)%NAME, BINVEN, NBD   , 
     *  QINDI,STRA41%R,STRA42%R,STRA43%R,STRA44%R,STSTOT%R,
     *  STSDER%R      , STOLD%R       , STNEW%R       , STRA31%R      ,
     *  STRA32%R      , STRA33%R      , STRA34%R      , STRA35%R      ,
     *  STRA36%R      , STRA37%R      , STRA38%R      , STRA39%R      ,
     *  ST1%R , ST2%R , ST3%R , ST4%R , STRA01%R      , SBETA%R  )
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE SEMIMP'
      ENDIF
C
C.....11.5 PASSAGE (EVENTUEL) EN FREQUENCE ABSOLUE.
C     """""""""""""""""""""""""""""""""""""""""""""
      IF (COURAN) THEN
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE TRANSF'
        CALL TRANSF
     *( STRA02%R      , SF%R  , SFR%R , SDFR%R, SCOSTE%R      ,
     *  SSINTE%R      , SUC%R , SVC%R , SXK%R , SITR11%I      ,
     *  SITR12%I      , SITR13%I      , STRA31%R      , STRA32%R      ,
     *  NPOIN2, NPLAN , NF    , RAISF , LT    , GRADEB, GRAPRD)
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE TRANSF'
      ELSE
        CALL OV ( 'X=Y     ' ,STRA02%R, SF%R , Z , C , NPOIN3*NF )
      ENDIF
C
C.....11.6 TEST POUR SAVOIR SI ON IMPRIME OU PAS.
C     """""""""""""""""""""""""""""""""""""""""""
      IMPRES=.FALSE.
      DEBRES=.FALSE.
      IF ((LT.GE.GRADEB).AND.(MOD(LT-GRADEB,GRAPRD).EQ.0)) IMPRES=.TRUE.
      IF (LT.EQ.GRADEB) DEBRES=.TRUE.
C
      IF (IMPRES) THEN
C
C.....11.7 IMPRESSION (EVENTUELLE) DES VARIABLES SUR LE MAILLAGE 2D.
C     """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
         IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE DUMP2D'
         CALL DUMP2D(LT, DEUPI, STRA02%R, NPTL )
         IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE DUMP2D'
C
         IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE BIEF_DESIMP'
         CALL BIEF_DESIMP(WAC_FILES(WACRES)%FMT,VARSOR,
     *            HIST,0,NPOIN2,WAC_FILES(WACRES)%LU,'STD',AT,
     *            LT,GRAPRD,GRAPRD,
     *            SORLEO,SORIMP,MAXVAR,TEXTE,GRADEB,GRADEB)
         IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE BIEF_DESIMP'
C
C.....11.8 IMPRESSION (EVENTUELLE) DES SPECTRES DIRECTIONNELS.
C     """"""""""""""""""""""""""""""""""""""""""""""""""""""""
         IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE ECRSPE'
         CALL ECRSPE
     *( STRA02%R, STRA01%R , STETA%R, NPLAN ,
     *  SFR%R   , NF  , NF   , NPOIN2      , AT ,LT ,
     *  STRA01%R, SITR01%I, NOLEO  , NPLEO, WAC_FILES(WACLEO)%LU ,
     *  WAC_FILES(WACLEO)%FMT , DEBRES, TITCAS, DATE  , TIME ,
     *  ISLEO , MESH%KNOLG%I )
         IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE ECRSPE'
C
      ENDIF
C
   10 CONTINUE
C
C=====C
C  12 C IMPRESSIONS GLOBALES (EVENTUELLES) EN FIN DE CALCUL.
C=====C=====================================================
C
      IF (GLOB)  THEN
        CALL IMPR(LISPRD,NIT,AT,NIT,6)
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE SOR3D'
        CALL SOR3D
     *( SF%R  , NPLAN  , NF       , STETA%R   , SFR%R ,
     *  NELEM2, NPOIN2 , AT       , SUC%R     , SVC%R ,
     *  SUV%R , SVV%R  , SDEPTH%R , VENT      , COURAN   ,
     *  MAREE , TITCAS , WAC_FILES(WACRBI)%LU , BINRBI    )
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE SOR3D'
      ENDIF
C
C----------------------------------------------------------------------
C
C  VALIDATION DES RESULTATS SUR LE FICHIER DE REFERENCES
C
      IF(VALID) CALL BIEF_VALIDA(BST1,TEXTE,
     *                   WAC_FILES(WACREF)%LU,WAC_FILES(WACREF)%FMT,
     *                   VARSOR,TEXTE,
     *                   WAC_FILES(WACRES)%LU,WAC_FILES(WACRES)%FMT,
     *                   MAXVAR,NPOIN2,NIT,NIT,ALIRE)
C
C----------------------------------------------------------------------
C
      RETURN
      END
