C                       ********************
                        SUBROUTINE TELEMAC2D
C                       ********************
C
     *(PASS,ATDEP,NITER,CODE,DTDEP,NEWTIME,DOPRINT)
C
C***********************************************************************
C  TELEMAC-2D VERSION 6.0 25/11/2009  J-M HERVOUET (LNHE) 01 30 87 80 18
C  VERSION AVEC TRACEURS MULTIPLES
C  DOUBLE APPEL DE SISYPHE (CHARRIAGE + SUSPENSION)
C  07/04/2008 : VARIATIONS EN TEMPS DES SOURCES
C  05/05/2008 : USIS ET VSIS DANS APPEL DE SISYPHE
C  14/05/2008 : HN INITIALISE POUR APPEL A SISYPHE
C  20/05/2008 : PASSAGE A TEL4DEL DES FLUX DUS AU LISSAGE DES VALEURS
C               NEGATIVES (NOUVEAU FILTER_H)
C  06/06/2008 : ARGUMENT OPTIONNEL BOUNDARY_COLOUR AJOUTE A LECLIM
C  16/06/2008 : SECOND APPEL A PROPIN PLACE APRES L'APPEL A BORD
C               (CHANGEMENTS DE CONDITIONS DANS BORD PAR L'UTILISATEUR)
C  25/06/2008 : DIFFIN2 REPREND SON NOM DIFFIN + ARGUMENT MESH
C  27/06/2008 : ARGUMENTS DE PROPIN_TELEMAC2D : MESH AJOUTE A LA FIN
C  29/07/2008 : APPEL A FLUSEC AJOUTE AVANT LE PREMIER CALL PRERES
C  13/08/2008 : APPEL A CHARAC CHANGE, ET CONDITIONS D'APPEL
C  20/08/2008 : LIST_PTS MODIFIED IN PARALLEL 
C  02/09/2008 : APPEL DE TEL4DEL MODIFIE (VITESSE ET DIFFUSION EN PLUS) 
C  25/09/2008 : APPEL DE TEL4DEL MODIFIE (FLUX ENVOYES PAR MESH%W%R) 
C  21/10/2008 : APPEL DE MASKTO MODIFIE (VERSION PARALLELE DE MASKTO)
C  09/02/2009 : SI CLIPPING DE H, AVEC HMIN AU LIEU DE 0.D0
C  16/02/2009 : APPEL DE POSITIVE_DEPTHS
C  19/02/2009 : IN CASE OF COMPUTATION CONTINUED, H CLIPPED
C  02/04/2009 : NEW FILE STRUCTURE T2D_FILES AND MED FORMAT
C  09/07/2009 : ARGUMENT NPTFR2 ADDED TO LECLIM
C  20/07/2009 : 1 OUT OF 3 CALLS TO TEL4DEL REMOVED (THANKS TO A
C               MODIFICATION OF PROPAG: COMPUTATION OF UDEL AND VDEL
C               IF(SOLSYS.EQ.1)
C  22/07/2009 : 3 NEW ARGUMENTS IN PROPAG
C***********************************************************************
C
C               TTTTT EEEEE L     EEEEE M   M   AA  CCCCC
C                 T   E     L     E     MM MM  A  A C
C                 T   EEE   L     EEE   M M M  AAAA C
C                 T   E     L     E     M   M  A  A C
C                 T   EEEEE LLLLL EEEEE M   M  A  A CCCCC
C
C         ******************************************************
C         *                                                    *
C         *            RESOLUTION DES EQUATIONS DE             *
C         *                   SAINT-VENANT                     *
C         *                EN VARIABLES U,V,H                  *
C         *                                                    *
C         ******************************************************
C
C
C      ADJO = .TRUE.  : MODE DIRECT
C      ADJO = .FALSE. : MODE ADJOINT
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |    PASS        | -->| -1 : ALL STEPS  
C |                |    |  0 : ONLY INITIALISATION 
C |                |    |  1 : ONLY TIME-STEPS STEPS  
C |    ATDEP       | -->| STARTING TIME WHEN CALLED FOR COUPLING
C |    NITER       | -->| NUMBER OF ITERATIONS WHEN CALLED FOR COUPLING
C |    CODE        | -->| CALLING PROGRAM (IF COUPLING) 
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C***********************************************************************
C
C SOUS-PROGRAMME APPELE PAR : HOMERE
C
C***********************************************************************
C
C-----------------------------------------------------------------------
C                    DECLARATION DES TYPES ET DIMENSIONS
C-----------------------------------------------------------------------
C
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_TELEMAC2D
      USE INTERFACE_TELEMAC2D, EX_TELEMAC2D => TELEMAC2D
      USE INTERFACE_SISYPHE, ONLY: SISYPHE
      USE GRACESTOP
      USE FRICTION_DEF
C     MODULE SPECIFIC TO COUPLING WITH ESTEL-3D
      USE M_COUPLING_ESTEL3D
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER,          INTENT(IN) :: PASS,NITER
      DOUBLE PRECISION, INTENT(IN) :: ATDEP
      CHARACTER(LEN=*), INTENT(IN) :: CODE
C     TIME STEP TO USE WHEN COUPLING WITH ESTEL-3D
      DOUBLE PRECISION, INTENT(IN), OPTIONAL :: DTDEP
C     ARE WE STARTING A NEW TIME STEP OR JUST ITERATING
      LOGICAL,          INTENT(IN), OPTIONAL :: NEWTIME
C     DO WE WANT TELEMAC2D TO OUTPUT IN THE LISTING OR NOT?
      LOGICAL,          INTENT(IN), OPTIONAL :: DOPRINT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C 
C ENTIERS
C
      INTEGER IELM,I,IELMX,ISOUSI,IBID,STOP2,LEOPRD_CHARR
      INTEGER ALIRE(MAXVAR),TROUVE(MAXVAR+10)
C
C SCALAIRES REELS
C
      DOUBLE PRECISION KMIN,KMAX,KARMAN,FLUSOR,FLUENT,HIST(1),AT0
      DOUBLE PRECISION C,MASSES,RELAXS,CFLMAX,TETAHC,DTCAS,RELAX
      DOUBLE PRECISION EMAX,EMIN,SCHMIT,ESTAR,SIGMAE,SIGMAK,C2,C1,CMU
C
C FOR TRACERS
C 
C     MASSOU: MASS CREATED BY SOURCE TERM DURING THE TIME STEP
C     MASTR0: INITIAL MASS
C     MASTR2: CURRENT MASS
C     MASTEN: MASS ENTERED THROUGH BOUNDARIES
C     MASTOU: TOTAL MASS CREATED BY SOURCE TERM    
      DOUBLE PRECISION MASSOU(MAXTRA),MASTR0(MAXTRA),MASTR2(MAXTRA)
      DOUBLE PRECISION MASTEN(MAXTRA),MASTOU(MAXTRA)
C
C LOGIQUES
C
      LOGICAL AKEP,INFOGS,INFOGT,ARRET1,ARRET2,YASMH,ARRET3,CORBOT
      LOGICAL CHARR,SUSP,SUSP1,NON,INIFLOW,YAFLODEL
      LOGICAL YASMI(MAXTRA)
C
      CHARACTER(LEN=24), PARAMETER :: CODE1='TELEMAC2D               '
      CHARACTER(LEN=16) :: FORMUL
C
C-----------------------------------------------------------------------
C
      INTEGER IOPTAN,IMAX,ITRAC,NPTFR2
C
C-----------------------------------------------------------------------
C
C AJOUTE POUR LES SCHEMAS CINETIQUES
C
      DOUBLE PRECISION FLUTSOR(MAXTRA),FLUTENT(MAXTRA),DTN
      DOUBLE PRECISION FLUSORTN,FLUENTN,TMAX,DTT
      INTEGER LTT
C
C-----------------------------------------------------------------------
C
C     FOR SISYPHE : GRAIN FEEDING AND CONSTANT FLOW DISCHARGE
      INTEGER :: ISIS_CFD, NSIS_CFD
      LOGICAL :: SISYPHE_CFD, CONSTFLOW_SIS
C     FRICTION DATA
      INTEGER :: KFROT_TP
C
      INTEGER  P_IMAX,P_IMIN
      DOUBLE PRECISION P_DMIN
      EXTERNAL P_IMAX,P_IMIN,P_DMIN
C
C     A ENLEVER SI ON MET USE INTERFACE_TELEMAC2D
!     DOUBLE PRECISION DEBSCE,TRSCE
!     EXTERNAL         DEBSCE,TRSCE
C
C-----------------------------------------------------------------------
C
      INTRINSIC MAX
C
C-----------------------------------------------------------------------
C
      DATA HIST /9999.D0/
C
C-----------------------------------------------------------------------
C
C  VARIABLES A LIRE EN CAS DE SUITE :
C  0 : ON NE LIT PAS    1 : ON LIT  (VOIR SS-PG NOMVAR)
C
C                                 0 : ANCIENNE PLACE DU TRACEUR
      DATA ALIRE /1,1,0,1,1,1,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     *            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     *            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     *            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
C
C-----------------------------------------------------------------------
C
C     CONVECTEUR POUR APPEL DE SISYPHE
C
      TYPE(BIEF_OBJ), POINTER :: USIS,VSIS
C
C-----------------------------------------------------------------------
C 
      SAVE
C
C-----------------------------------------------------------------------
C
      NON=.FALSE.
      CHARR=.FALSE.
      SUSP=.FALSE.
C
C-----------------------------------------------------------------------
C
C     FOR INITIALISATION OF FLODEL (ARRAY FLOW) IN DELWAQ
C
      IF(OPTBAN.EQ.1.OR.OPTBAN.EQ.3) THEN
        INIFLOW=.FALSE.
      ELSE
        INIFLOW=.TRUE.
      ENDIF
C     
C     FOR COMPUTING EXTRA FLOWS DUE TO TIDAL FLATS TREATMENT
C
      IF(INCLUS(COUPLING,'DELWAQ')) THEN
        YAFLODEL=.TRUE.
      ELSE
        YAFLODEL=.FALSE.
      ENDIF
C
      IF(NTRAC.GT.0) THEN
        DO ITRAC=1,NTRAC
C         VOIR POINT_TELEMAC2D
          ALIRE(31+ITRAC) = 1
        ENDDO
      ENDIF
C
      IF(ITURB.NE.3) ALIRE(10) = 0
      IF(ITURB.NE.3) ALIRE(11) = 0
      IF(ITURB.EQ.1) ALIRE(12) = 0
C
C-----------------------------------------------------------------------
C
C     USE DOPRINT TO LIMIT TELEMAC-2D OUTPUTS IN THE LISTING
      IF(PRESENT(DOPRINT)) THEN
        LISTIN =  DOPRINT
        ENTET  =  DOPRINT
      ENDIF
C
C-----------------------------------------------------------------------
C
      IF(PASS.EQ.0) THEN
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'INITIALISATION DE TELEMAC2D POUR ',CODE
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) 'INITIALISING TELEMAC2D FOR ',CODE
        ENDIF
      ELSEIF(PASS.EQ.1) THEN
        GO TO 700
      ELSEIF(PASS.NE.-1) THEN
        IF(LNG.EQ.1) WRITE(LU,*) 'MAUVAIS ARGUMENT PASS : ',PASS
        IF(LNG.EQ.2) WRITE(LU,*) 'WRONG ARGUMENT PASS: ',PASS
        CALL PLANTE(1)
        STOP
      ENDIF
C
C=======================================================================
C
C : 1          LECTURE , PREPARATION ET CONTROLE DES DONNEES
C
C=======================================================================
C
C  TYPES DE DISCRETISATION : POUR L'INSTANT TRIANGLES P1
C
      IELM=IELM1
C     ELEMENT LE PLUS COMPLEXE
      IELMX = MAX(IELMH,IELMU,IELMT,IELMK,IELME)
C
C-----------------------------------------------------------------------
C
C LECTURE DES CONDITIONS LIMITES ET INDICES DES POINTS FRONTIERES.
C
      IF(IELMX.EQ.13) THEN
        NPTFR2=2*NPTFR
      ELSE
        NPTFR2=NPTFR
      ENDIF
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE LECLIM'
      CALL LECLIM (LIHBOR%I   , LIUBOR%I , LIVBOR%I , LITBOR%ADR(1)%P%I,
     *             HBOR%R     , UBOR%R   , VBOR%R   , TBOR%ADR(1)%P%R ,
     *             AUBOR%R    , ATBOR%ADR(1)%P%R   , BTBOR%ADR(1)%P%R ,
     *             MESH%NPTFR , 3        ,NTRAC.GT.0,
     *             T2D_FILES(T2DCLI)%LU,
     *             KENT       , KENTU    , KSORT ,  KADH , KLOG , KINC,
     *             NUMLIQ%I   ,MESH,BOUNDARY_COLOUR%I,NPTFR2)
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE LECLIM'
C
C DUPLICATION DES CONDITIONS AUX LIMITES DES TRACEURS
C
      IF(NTRAC.GE.2) THEN
        DO ITRAC=2,NTRAC
          DO I=1,NPTFR
            LITBOR%ADR(ITRAC)%P%I(I)=LITBOR%ADR(1)%P%I(I)
              TBOR%ADR(ITRAC)%P%R(I)=  TBOR%ADR(1)%P%R(I)
             ATBOR%ADR(ITRAC)%P%R(I)= ATBOR%ADR(1)%P%R(I)
             BTBOR%ADR(ITRAC)%P%R(I)= BTBOR%ADR(1)%P%R(I)
          ENDDO
        ENDDO
      ENDIF
C
C-----------------------------------------------------------------------
C  COMPLEMENT DE LA STRUCTURE DE DONNEES POUR BIEF
C-----------------------------------------------------------------------
C
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE INBIEF'
      CALL INBIEF(LIHBOR%I,KLOG,IT1,IT2,IT3,LVMAC,IELMX,
     *            LAMBD0,SPHERI,MESH,T1,T2,OPTASS,PRODUC,EQUA)
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE INBIEF'
C
      IF(IELMX.EQ.13) THEN
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE COMPLIM'
        CALL COMPLIM( LIUBOR%I , LIVBOR%I , LITBOR%ADR(1)%P%I,
     *                UBOR%R   , VBOR%R   , TBOR%ADR(1)%P%R ,
     *                AUBOR%R    , ATBOR%ADR(1)%P%R , BTBOR%ADR(1)%P%R ,
     *                MESH%NBOR%I,MESH%NPTFR , MESH%NPOIN, NTRAC.GT.0, 
     *                KENT , KENTU , KSORT ,KADH , KLOG , KINC,
     *                IELMU,IELMU,IELMT,MESH)
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE COMPLIM'
      ENDIF 
C
C-----------------------------------------------------------------------
C  DEFINITION DE ZONES PAR L'UTILISATEUR
C-----------------------------------------------------------------------
C
      IF(DEFZON) CALL DEF_ZONES
C
C-----------------------------------------------------------------------
C  CHANGING GLOBAL TO LOCAL IN LIST OF POINTS IN PARALLEL
C-----------------------------------------------------------------------
C
      IF(NPTS.GT.0.AND.NCSIZE.GT.0) THEN
        DO I=1,NPTS
          LIST_PTS(I)=MESH%KNOGL%I(LIST_PTS(I))
        ENDDO
      ENDIF
C
C-----------------------------------------------------------------------
C  RECHERCHE DU FOND ET DU FROTTEMENT DANS LE FICHIER DE GEOMETRIE:
C-----------------------------------------------------------------------
C
      IF(     .NOT.INCLU2(ESTIME,'FROTTEMENT')
     *   .AND..NOT.INCLU2(ESTIME,'FRICTION'  )  ) THEN
C       NO PARAMETER ESTIMATION
        CALL FONSTR(T1,ZF,T2,CHESTR,T2D_FILES(T2DGEO)%LU,
     *              T2D_FILES(T2DFON)%LU,T2D_FILES(T2DFON)%NAME,
     *              MESH,FFON,LISTIN)
        CORBOT=.TRUE.
      ELSEIF(NITERA.EQ.1.AND..NOT.ADJO) THEN
C       WITH PARAMETER ESTIMATION (HENCE NITERA DEFINED), 
C       FONSTR CALLED ONCE TO GET
C       THE BOTTOM TOPOGRAPHY AND THE INITIAL FRICTION (CALL TO STRCHE)
        CALL FONSTR(T1,ZF,T2,CHESTR,T2D_FILES(T2DGEO)%LU,
     *              T2D_FILES(T2DFON)%LU,T2D_FILES(T2DFON)%NAME,
     *              MESH,FFON,LISTIN)
C       IF OPTID=0, VALUES OF SETSTR GIVEN BY FILE, MUST NOT BE ERASED
        IF(OPTID.NE.0) CALL INITSTR(CHESTR,SETSTR,ZONE%I,NZONE,NPOIN,T1)
        CALL ASSIGNSTR(CHESTR,SETSTR,ZONE%I,NZONE,NPOIN)
        CORBOT=.TRUE.
      ELSE
C       IN PARAMETER ESTIMATION, FROM NITERA=2 ON, BOTTOM IS NOT READ
C       AGAIN, SO NO CALL TO CORFON
        CORBOT=.FALSE.     
      ENDIF
C
C     INITIALIZATION OF FRICTION COEFFICIENT BY ZONE
C
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE FRICTION_CHOICE' 
      CALL FRICTION_CHOICE(0,KARMAN)
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE FRICTION_CHOICE'
C
C-----------------------------------------------------------------------
C
C PREPARATION DU FICHIER DE RESULTATS (FACULTATIF)
C
C      STANDARD SELAFIN
C
      IF(ADJO) THEN
C
        IF(T2D_FILES(T2DRBI)%NAME.NE.' '.AND.
     *     INCLU2(ESTIME,'DEBUG')) THEN
         CALL ECRGEO(MESH%X%R,MESH%Y%R,MESH%NPOIN,MESH%NBOR%I,
     *               T2D_FILES(T2DRBI)%LU,IBID,TEXTE,VARCLA,NVARCL,
     *               TITCAS,SORLEOA,MAXVAR,MESH%IKLE%I,
     *               MESH%NELEM,MESH%NPTFR,3,MARDAT,MARTIM,
     *               NCSIZE,NPTIR,MESH%KNOLG%I,I3=I_ORIG,I4=J_ORIG)
        ENDIF
C
      ELSE
C
!       CALL ECRGEO(MESH%X%R,MESH%Y%R,MESH%NPOIN,MESH%NBOR%I,
!    *            NRES,IBID,TEXTE,VARCLA,NVARCL,
!    *            TITCAS,SORLEO,MAXVAR,MESH%IKLE%I,
!    *            MESH%NELEM,MESH%NPTFR,3,MARDAT,MARTIM,
!    *            NCSIZE,NPTIR,MESH%KNOLG%I,I3=I_ORIG,I4=J_ORIG)
        ! CREATION DU JEU DE DONNEES POUR UN FORMAT DE FICHIER
        ! FORMAT_RES.
        ! LE JEU DE DONNEES EST CREE DANS LE FICHIER NRES, ET EST
        ! DEFINIT PAR UN TITRE ET LES VARIABLES A ECRIRE. 
        CALL CREATE_DATASET(T2D_FILES(T2DRES)%FMT, ! FORMAT FICHIER RESULTAT
     *                      T2D_FILES(T2DRES)%LU,  ! LU FICHIER RESULTAT
     *                      TITCAS,     ! TITRE DE L'ETUDE
     *                      MAXVAR,     ! MAX VARIABLES SORTIE
     *                      TEXTE,      ! NOMS VARIABLES SORTIE
     *                      SORLEO)     ! SORTIE OU PAS DES VARIABLES
        ! ECRITURE DU MAILLAGE DANS LE FICHIER SORTIE :
        ! SI ON EST ON PARALLEL, FAUT L'INDIQUER VIA NCSIZE ET NPTIR.
        ! LES AUTRES INFORMATIONS SONT DANS MESH.
        ! EN PLUS : DATE/TEMPS DE DEPART ET LES COORDONNEES DE
        ! L'ORIGINE.
        CALL WRITE_MESH(T2D_FILES(T2DRES)%FMT, ! FORMAT FICHIER RESULTAT     
     *                  T2D_FILES(T2DRES)%LU,  ! LU FICHIER RESULTAT
     *                  MESH,          ! DESCRIPTEUR MAILLAGE
     *                  1,             ! NOMBRE DE PLAN /NA/
     *                  MARDAT,        ! DATE DEBUT
     *                  MARTIM,        ! HEURE DEBUT
     *                  I_ORIG,J_ORIG) ! COORDONNEES DE L'ORIGINE.
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C     INITIALISATION DU BLOC PRIVE
C
      IF(NPRIV.GT.0) CALL OS('X=0     ',X=PRIVE)
C
C  ON RALLONGE COSLAT ET SINLAT POUR LEUR DONNER LA DIMENSION DE U ET V
C  IDEM POUR LE FROTTEMENT
C
      IF(IELMU.NE.IELM1) THEN
        IF(SPHERI) CALL CHGDIS(MESH%COSLAT,IELM1,IELMU,MESH)
        IF(SPHERI) CALL CHGDIS(MESH%SINLAT,IELM1,IELMU,MESH)
        CALL CHGDIS(CHESTR,IELM1,IELMU,MESH)
      ENDIF
C
C=======================================================================
C
C  REPERAGE DES FRONTIERES
C
      IF(NCSIZE.GT.1) THEN
       NFRLIQ=0
       DO I=1,NPTFR
         NFRLIQ=MAX(NFRLIQ,NUMLIQ%I(I))
       ENDDO
       NFRLIQ=P_IMAX(NFRLIQ)
       WRITE(LU,*) ' '
       IF(LNG.EQ.1) WRITE(LU,*) 'NOMBRE DE FRONTIERES LIQUIDES :',NFRLIQ
       IF(LNG.EQ.2) WRITE(LU,*) 'NUMBER OF LIQUID BOUNDARIES:',NFRLIQ
      ELSE
       IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE FRONT2'
       CALL FRONT2(NFRLIQ,NFRSOL,DEBLIQ,FINLIQ,DEBSOL,FINSOL,
     *             LIHBOR%I,LIUBOR%I,
     *             MESH%X%R,MESH%Y%R,MESH%NBOR%I,MESH%KP1BOR%I,
     *             IT1%I,NPOIN,NPTFR,KLOG,LISTIN,NUMLIQ%I,MAXFRO)
       IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE FRONT2'
      ENDIF
C
C=======================================================================
C
C  LECTURE DU FICHIER DES COURBES DE TARAGE
C
      IF(T2D_FILES(T2DMAB)%NAME(1:1).NE.' ') THEN
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE READ_FIC_CURVES'
        CALL READ_FIC_CURVES(T2D_FILES(T2DMAB)%LU,NFRLIQ,
     *                       STA_DIS_CURVES,PTS_CURVES)
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE READ_FIC_CURVES'
      ENDIF
C
C=======================================================================
C
C CORRECTION DES NORMALES AUX NOEUDS DE BORD POUR AVOIR DES NORMALES
C AU SEGMENT LIQUIDE ADJACENT DANS LE CAS DES TRANSITIONS LIQUIDE-SOLIDE
C
      IF(EQUA(1:15).NE.'SAINT-VENANT VF') THEN
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE CORNOR'
        CALL CORNOR(MESH%XNEBOR%R,MESH%YNEBOR%R,
     *              MESH%XSGBOR%R,MESH%YSGBOR%R,
     *              MESH%KP1BOR%I,NPTFR,KLOG,LIHBOR%I,
     *              T1,T2,MESH)
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE CORNOR'
      ENDIF
C
C=======================================================================
C
C REMPLISSAGE PAR DEFAUT DU TABLEAU MASKEL
C (TOUS LES ELEMENTS SONT A CONSIDERER)
C
      IF(MSK) CALL OS ( 'X=C     ' , MASKEL , S , S , 1.D0 )
C
C     MASQUAGE ARTIFICIEL DONNE PAR L'UTILISATEUR
C     CE SOUS-PROGRAMME EST AUSSI APPELE A CHAQUE PAS DE TEMPS
      IF(MSKUSE) THEN
        CALL MASKOB (MASKEL%R,MESH%X%R,MESH%Y%R,
     *               IKLE%I,NELEM,NELMAX,NPOIN,0.D0,0)
      ENDIF
C
C-----------------------------------------------------------------------
C  INTEGRAL OF TEST FUNCTIONS (ONCE FOR ALL AND WITHOUT MASKING) 
C-----------------------------------------------------------------------
C
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE MASBAS2D'
      CALL MASBAS2D(VOLU2D,V2DPAR,UNSV2D,IELM1,MESH,.FALSE.,
     *              MASKEL,T2,T2)
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE MASBAS2D'
C
C=======================================================================
C
C CORRECTING THE BOTTOM WITH USER-SUBROUTINE CORFON
C ZF IS TREATED AS LINEAR IN CORFON
C IF(CORBOT) : SEE CALL FONSTR ABOVE, IN PARAMETER ESTIMATION, ZF
C IS READ ONLY AT THE FIRST RUN
C
      IF(CORBOT) THEN 
        IF(IELMH.NE.IELM1) CALL CHGDIS(ZF,IELMH,IELM1,MESH)
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE CORFON'
        CALL CORFON
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE CORFON'
        IF(IELMH.NE.IELM1) CALL CHGDIS(ZF,IELM1,IELMH,MESH)
      ENDIF
C
C=======================================================================
C
C REDEFINITION EVENTUELLE DES CARACTERISTIQUES DES SOURCES
C
C EN STANDARD, SOURCE NE FAIT RIEN
C
      CALL SOURCE_TELEMAC2D
C
C=======================================================================
C
C ANALYSE FINE DE LA TOPOGRAPHIE
C
      IF(OPTBAN.EQ.2) THEN
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE TOPOGR'
        CALL TOPOGR(ZF%R,T1%R,ZFE%R,IKLE%I,MESH%IFABOR%I,
     *              MESH%NBOR%I,MESH%NELBOR%I,MESH%NULONE%I,
     *              IT1%I,IT2%I,IT3%I,
     *              NELEM,NPTFR,NPOIN,MXPTVS)
       IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE TOPOGR'
      ENDIF
C
C=======================================================================
C
C : 2                  INITIALISATIONS
C
C=======================================================================
C
C CONSTANTES K-EPSILON ET AUTRES (KARMAN EST UTILISEE MEME SANS K-E)
C
      CALL COSAKE(KARMAN,CMU,C1,C2,SIGMAK,SIGMAE,
     *            ESTAR,SCHMIT,KMIN,KMAX,EMIN,EMAX)
C
      IF(ITURB.EQ.3) THEN
C       IL FAUDRA INITIALISER K ET EPSILON
        AKEP = .TRUE.
      ELSE
C       IL NE FAUDRA PAS INITIALISER K ET EPSILON
        AKEP = .FALSE.
      ENDIF
C
C INITIALISATION DES GRANDEURS PHYSIQUES
C
C     CONDIN EST APPELE MEME EN CAS DE SUITE, POUR AVOIR UNE DEFINITION
C     DE C0 QUI NE CHANGE PAS EN CAS DE SUITE (CAS DES ONDES INCIDENTES)
C
      IF(ADJO) THEN
        CALL CONDIN_ADJ(ALIRE,T2D_FILES(T2DRES)%LU,TROUVE)
      ELSE
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE CONDIN'
        CALL CONDIN
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE CONDIN'
      ENDIF
C
C     ON CORRIGE LES ERREURS DE L'UTILISATEUR SI IL A MIS H<0
C     ICI NOMBRE DE POINTS FORCE A NPOIN.
      CALL CLIP(H,0.D0,.TRUE.,1.D6,.FALSE.,-NPOIN)
C
C     STOCKAGE DE LA CELERITE INITIALE (POUR L'ONDE INCIDENTE)
C
      CALL CELERITE
C
C CALCUL DE LA HAUTEUR DE REFERENCE POUR LES EQUATIONS DE BOUSSINESQ
C
      IF(EQUA(1:10).EQ.'BOUSSINESQ') THEN
        CALL HREF
      ENDIF
C     AJOUTE LE 27/05/2002 (AVANT : NON INITIALISE)
      AT0=0.D0
C
      IF(.NOT.DEBU.AND..NOT.ADJO) THEN
C
C       BEWARE : BIEF_SUITE WILL TAKE THE BOTTOM IN THE FILE
C                IF IT IS THERE.
C
C       FRICTION COEFFICIENT ALSO READ IN CASE IT HAS BEEN DONE
C       BY THE USER INTERFACE (JMH 27/11/2006)
        ALIRE(19)=1
        CALL BIEF_SUITE(VARSOR,VARCL,IBID,
     *                  T2D_FILES(T2DPRE)%LU,
     *                  T2D_FILES(T2DPRE)%FMT,
     *                  HIST,0,NPOIN,AT,TEXTPR,VARCLA,
     *                  NVARCL,TROUVE,ALIRE,LISTIN,.TRUE.,MAXVAR)
        ALIRE(19)=0
        IF(RAZTIM) THEN
          AT=0.D0
          IF(LNG.EQ.1) WRITE(LU,*) 'TEMPS ECOULE REMIS A ZERO'
          IF(LNG.EQ.2) WRITE(LU,*) 'ELAPSED TIME RESET TO ZERO'
        ENDIF
        AT0=AT
        CALL RESCUE(U%R,V%R,H%R,FV%R,ZF%R,T,TRAC0,NTRAC,
     *              ITURB,NPOIN,AKEP,TROUVE)
C       CAS OU HAUTEURS POSITIVES NECESSAIRE
        IF(OPTBAN.EQ.1.AND.OPT_HNEG.EQ.2) THEN
          CALL CLIP(H,0.D0,.TRUE.,1.D6,.FALSE.,-NPOIN)
        ENDIF
      ENDIF
C
      TMAX=DUREE+AT0
C
C-----------------------------------------------------------------------
C
C  INITIALISATIONS PROPRES AU VOLUMES FINIS
C
C-----------------------------------------------------------------------
C
      IF(EQUA(1:15).EQ.'SAINT-VENANT VF') THEN
C
        CALL OS( 'X=YZ    ' , QU , U , H , C )
        CALL OS( 'X=YZ    ' , QV , V , H , C )
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      LT=0
      LTT=0
C
C=======================================================================
C EXTENSION DES VARIABLES QUI NE SONT PAS LINEAIRES P1
C=======================================================================
C
      IF(NTRAC.GT.0.AND.IELMT.NE.IELM1) THEN
        DO ITRAC=1,NTRAC
          CALL CHGDIS( T%ADR(ITRAC)%P ,IELM1 , IELMT , MESH )
        ENDDO
      ENDIF
      IF(IELMH.NE.IELM1) THEN
        CALL CHGDIS( H  , IELM1 , IELMH , MESH )
        CALL CHGDIS( ZF , IELM1 , IELMH , MESH )
      ENDIF
      IF(IELMU.NE.IELM1) THEN
        CALL CHGDIS( U , IELM1 , IELMU , MESH )
        CALL CHGDIS( V , IELM1 , IELMU , MESH )
      ENDIF
C
C=======================================================================
C CONDITIONS INITIALES NON CONTENUES DANS SUITE OU CONDIN
C=======================================================================
C
C  CLIPPING (CONDITIONNEL) DE H
C
      IF(CLIPH) CALL CLIP( H , HMIN , .TRUE. , 1.D6 , .FALSE. , 0 )
C
C-----------------------------------------------------------------------
C CONDITIONS INITIALES DE METEO
C
      IF (VENT.OR.ATMOS) THEN
        CALL METEO(PATMOS%R,WINDX%R,WINDY%R,
     *             FUAIR,FVAIR,MESH%X%R,MESH%Y%R,AT,LT,NPOIN,VENT,ATMOS,
     *             H%R,T1%R,GRAV,ROEAU,NORD,PRIVE)
      ENDIF
C
C-----------------------------------------------------------------------
C
C LECTURE DE LA GEOMETRIE DES SINGULARITES
C
      IF(NWEIRS.GT.0) THEN
       CALL LECSNG(NWEIRS,NWRMAX,NPSING,NUMDIG%I,
     *             ZDIG%R,PHIDIG%R,IOPTAN,NPSMAX,NPOIN,
     *             T2D_FILES(T2DFO1)%LU)
      ENDIF
      IF(NSIPH.GT.0) THEN
       CALL LECSIP(RELAXS,NSIPH,ENTSIP,SORSIP,SECSCE,
     *             ALTSCE,CSSCE,CESCE,DELSCE,
     *             ANGSCE,LSCE,MAXSCE,T2D_FILES(T2DFO1)%LU)
      ENDIF
C
C-----------------------------------------------------------------------
C
C CONDITIONS INITIALES POUR LE MODELE K-EPSILON ET LA DIFFUSION
C
C   K-EPSILON
C
C     SI AKEP = .FALSE. K ET EPSILON ONT ETE FAITS DANS SUITE OU CONDIN
      IF(AKEP) THEN
C
        CALL FRICTION_CHOICE(1, KARMAN)
        IF(FRICTB) THEN
           KFROT_TP = 0
           IF(KFROT.EQ.NZONES) KFROT_TP = 1 ! NEED A NON ZERO VALUE
        ELSE
           KFROT_TP = KFROT
        ENDIF
!
        CALL AKEPIN(AK%R,EP%R,U%R,V%R,H%R,NPOIN,KFROT_TP,CMU,C2,
     *              ESTAR,SCHMIT,KMIN,EMIN,CF%R)
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C     PREPARATION DES CONDITIONS AUX LIMITES POUR LES SEUILS.
C
      IF(NWEIRS.GT.0) THEN
C
        CALL CLSING(NWEIRS,NPSING,NPSMAX,NUMDIG%I,
     *              MESH%X%R,MESH%Y%R,ZF%R,CHESTR%R,NKFROT%I,
     *              KARMAN,ZDIG%R,PHIDIG%R,MESH%NBOR%I,
     *              H%R,T,NTRAC,IOPTAN,T1%R,UBOR%R,VBOR%R,TBOR,
     *              LIHBOR%I,LIUBOR%I,LIVBOR%I,LITBOR,GRAV)
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C     TYPES DE CONDITIONS POUR LE TRACEUR :
C
      IF(NTRAC.GT.0) THEN
        IF(NWEIRS.GT.0) CALL CLTRAC(NWEIRS,NPSING,NPSMAX,NUMDIG%I,
     *                 ZF%R,ZDIG%R,H%R,T,MESH%NBOR%I,LITBOR,TBOR,NTRAC)
        DO ITRAC=1,NTRAC
        CALL DIFFIN(MASKTR,LIMTRA%I,LITBOR%ADR(ITRAC)%P%I,
     *              IT1%I,U%R,V%R,MESH%XNEBOR%R,MESH%YNEBOR%R,
     *              MESH%NBOR%I,MESH%KP1BOR%I,NPTFR,
     *              KENT,KSORT,KLOG,KINC,KNEU,KDIR,KDDL,
     *              ICONVF(3),MESH%NELBOR%I,NPOIN,NELMAX,MSK,MASKEL%R,
     *              NFRLIQ,THOMFR,DEBLIQ,FINLIQ,FRTYPE,
     *              TN%ADR(ITRAC)%P,TBOR%ADR(ITRAC)%P,MESH)
        ENDDO
      ENDIF
C
C     TYPES DE CONDITIONS POUR LA PROPAGATION:
C
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE PROPIN'
      CALL PROPIN_TELEMAC2D
     *            (LIMPRO%I,LIMPRO%DIM1,MASK,LIUBOR%I,LIVBOR%I,
     *             LIHBOR%I,MESH%KP1BOR%I,MESH%NBOR%I,NPTFR,
     *             KENT,KENTU,KSORT,KADH,KLOG,KINC,
     *             KNEU,KDIR,KDDL,KOND,CLH%I,CLU%I,CLV%I,
     *             U%ELM,U%R,V%R,GRAV,H%R,LT,NPOIN,
     *             MESH%NELBOR%I,NELMAX,MSK,MASKEL%R,
     *             NFRLIQ,THOMFR,DEBLIQ,FINLIQ,FRTYPE,
     *             MESH%XNEBOR%R,MESH%YNEBOR%R,ENTET,MESH)
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE PROPIN'
C
C     PROPIN SERA RAPPELE DANS LA BOUCLE EN TEMPS APRES CHAQUE APPEL
C     DE BORD
C
C-----------------------------------------------------------------------
C
C     COEFFICIENT DE FROTTEMENT :
C
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE FRICTION_CHOICE'
      CALL FRICTION_CHOICE(1,KARMAN)
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE FRICTION_CHOICE'
C
C  DIFFUSION DES VITESSES (APPELE ICI POUR INITIALISER VISC AU CAS
C                          OU C'EST UNE VARIABLE DEMANDEE EN SORTIE)
      IF(ITURB.EQ.1) THEN
C
        CALL OS('X=C     ', X=VISC , C=PROPNU )
C
      ELSEIF(ITURB.EQ.2) THEN
C
        CALL DISPER( VISC , U%R , V%R , H%R , CF%R , ELDER , PROPNU )
C
      ELSEIF(ITURB.EQ.3) THEN
C
        CALL VISTUR(VISC,AK,EP,NPOIN,CMU,PROPNU)
C
      ELSEIF(ITURB.EQ.4) THEN
C
        CALL SMAGOR(VISC,CF,U,V,MESH,T1,T2,T3,T4,MSK,MASKEL,PROPNU)
C
      ELSE
        IF(LISTIN) THEN
          IF(LNG.EQ.1) WRITE(LU,15) ITURB
          IF(LNG.EQ.2) WRITE(LU,16) ITURB
        ENDIF
15      FORMAT(1X,'ITURB=',1I6,'MODELE DE TURBULENCE NON PREVU')
16      FORMAT(1X,'ITURB=',1I6,'UNKNOWN TURBULENCE MODEL')
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C  FLOTTEUR(S)
C
      IF(NFLOT.NE.0) CALL FLOT(XFLOT%R,YFLOT%R,NFLOT,NITFLO,
     *                         FLOPRD,MESH%X%R,MESH%Y%R,
     *                         NPOIN,DEBFLO%I,FINFLO%I,NIT)
C
C-----------------------------------------------------------------------
C  DERIVE(S) LAGRANGIENNE(S)
C
      IF(NLAG.NE.0) CALL LAGRAN(NLAG,DEBLAG%I,FINLAG%I)
C
C-----------------------------------------------------------------------
C  POSITIONS DES POINTS DE REJET
C
      IF(NREJET.NE.0.OR.NREJTR.NE.0) THEN
        CALL PROXIM(ISCE,XSCE,YSCE,
     *              MESH%X%R,MESH%Y%R,
     *              MAX(NREJET,NREJTR),NPOIN,
     *              MESH%IKLE%I,NELEM,NELMAX)
      ENDIF
C
C=======================================================================
C FIN DES CONDITIONS INITIALES
C=======================================================================
C
C INITIALISATION DU MODULE DE CONVECTION
C FTILD CONTIENT UTILD,VTILD,HTILD,(TTILD),(AKTILD ET EPTILD)
C
      CALL OS( 'X=0     ' , X=FTILD )
C
C***********************************************************************
C
C IMPRESSIONS ET SORTIES POUR LES CONDITIONS INITIALES.
C
      IF(LISTIN) CALL ENTETE(1,AT,LT)
C
C     OUTINI IS KEY-WORD "OUTPUT OF INITIAL CONDITONS"
C     IT HAS PRIORITY OVER FIRST TIME-STEP FOR GRAPHIC PRINTOUTS.
C
C     NOTE THAT OUTPUTS ARE DONE WITHIN ESTEL3D IN COUPLED MODE)
C
      IF(OUTINI .AND. (.NOT.ADJO)
     *          .AND. (CODE(1:7).NE.'ESTEL3D') ) THEN
C
C SECTIONS DE CONTROLE (AVEC 0. A LA PLACE DE DT)
C
        IF(NCP.NE.0.AND.(ENTET.OR.CUMFLO)) THEN 
          CALL FLUSEC_TELEMAC2D(U,V,H,MESH%IKLE%I,MESH%XEL%R,MESH%YEL%R,
     *                          MESH%NELMAX,MESH%NELEM,
     *                          MESH%X%R,MESH%Y%R,
     *                          0.D0,NCP,CTRLSC,ENTET,AT,MESH%KNOGL%I,
     *                          MSKSEC,BM1,BM2,T1,H,MESH,S,CV1,
     *                          MESH%IFABOR%I,COMFLU,CUMFLO)
        ENDIF
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE PRERES_TELEMAC2D'
        CALL PRERES_TELEMAC2D
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE PRERES_TELEMAC2D'
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE DESIMP'
        CALL BIEF_DESIMP(T2D_FILES(T2DRES)%FMT,VARSOR,
     *                  HIST,0,NPOIN,T2D_FILES(T2DRES)%LU,'STD',AT,LT,
     *                  LISPRD,LEOPRD,
     *                  SORLEO,SORIMP,MAXVAR,TEXTE,0,     0)
C                                                  PTINIG,PTINIL
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE DESIMP'
C
      ENDIF
C
C=======================================================================
C
C     COUPLING WITH DELWAQ
C
      IF(INCLUS(COUPLING,'DELWAQ')) THEN
C
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE TEL4DEL'
C
C     T3 : MODIFIED DEPTH TO TAKE INTO ACCOUNT MASS-LUMPING
C          IN THE CONTINUITY EQUATION
      IF(ABS(1.D0-AGGLOC).GT.1.D-8) THEN
        CALL VECTOR(T3 ,'=','MASVEC          ',IELMH,
     *              1.D0-AGGLOC,H ,S,S,S,S,S,MESH,MSK,MASKEL)
        IF(NCSIZE.GT.1) CALL PARCOM(T3,2,MESH)
        CALL OS('X=XY    ',X=T3 ,Y=UNSV2D)
        CALL OS('X=X+CY  ',X=T3 ,Y=H ,C=AGGLOC)
      ELSE
        CALL OS('X=Y     ',X=T3 ,Y=H )
      ENDIF
      CALL TEL4DEL(MESH%NPOIN,
     *        MESH%NPOIN,MESH%NELEM,MESH%NSEG,MESH%IKLE%I,MESH%ELTSEG%I,
     *        MESH%GLOSEG%I,MESH%ORISEG%I,MESH%GLOSEG%DIM1,
     *        MESH%X%R,MESH%Y%R,MESH%NPTFR,LIHBOR%I,
     *        MESH%NBOR%I,1,AT,DT,LT,NIT,T3%R,H%R,T3%R,U%R,V%R,
     *        T%ADR(MAX(IND_S,1))%P%R,
     *        T%ADR(MAX(IND_T,1))%P%R,VISC%R,TITCAS,
     *        T2D_FILES(T2DGEO)%NAME,T2D_FILES(T2DCLI)%NAME,WAQPRD,
     * T2DDL1,T2D_FILES(T2DDL1)%NAME,T2DDL2,T2D_FILES(T2DDL2)%NAME,
     * T2DDL3,T2D_FILES(T2DDL3)%NAME,T2DDL5,T2D_FILES(T2DDL5)%NAME,
     * T2DDL6,T2D_FILES(T2DDL6)%NAME,T2DDL7,T2D_FILES(T2DDL7)%NAME,
     * T2DL11,T2D_FILES(T2DL11)%NAME,T2DDL4,T2D_FILES(T2DDL4)%NAME,
     * T2DDL8,T2D_FILES(T2DDL8)%NAME,T2DDL9,T2D_FILES(T2DDL9)%NAME,
     * T2DL10,T2D_FILES(T2DL10)%NAME,INFOGR,NELEM,SALI_DEL,TEMP_DEL,
     * VELO_DEL,DIFF_DEL,MARDAT,MARTIM,FLODEL%R,INIFLOW,MESH%W%R,
     * .FALSE.,FLULIM%R,V2DPAR%R,MESH%KNOLG%I,MESH,MESH)
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE TEL4DEL'
C
      ENDIF
C
C=======================================================================
C
C     OPTIONAL USER OUTPUT (COURTESY JACEK JANKOWSKI, BAW)
      CALL UTIMP_TELEMAC2D(LT,AT,PTINIG,LEOPRD,PTINIL,LISPRD)
C
C=======================================================================
C
C  INITIALISATION DES CHAMPS CONVECTEURS ET PROPAGATEURS
C
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE HPROPA'
      CALL HPROPA(HPROP,H,H,PROLIN,HAULIN,TETAC,NSOUSI)
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE HPROPA APPEL DE CHPCON'
      CALL CHPCON(UCONV,VCONV,U,V,U,V,TETAU)
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE CHPCON'
      IF(SOLSYS.EQ.2) THEN
C       INITIALISATION OF UDEL AND VDEL ONLY TO FIRST CALL TO SISYPHE
        CALL OS('X=Y     ',X=UDEL,Y=UCONV)
        CALL OS('X=Y     ',X=VDEL,Y=VCONV)
        USIS=>UDEL
        VSIS=>VDEL
      ELSE
        USIS=>UCONV
        VSIS=>VCONV
      ENDIF      
C
C=======================================================================
C
C     TETAHC : SEMI-IMPLICITATION DE H DANS L'EQUATION DE CONTINUITE
C              EST AUSSI UTILISE POUR LES FLUX DANS LE BILAN DE MASSE
      TETAHC = TETAC
      IF(ICONVF(2).EQ.5) TETAHC = 0.D0
C
C     FIRST COMPUTATION OF POROSITY
C
      IF(OPTBAN.EQ.3) THEN
        CALL POROS(TE5,ZF,H,MESH)
        IF(MSK) CALL OS('X=XY    ',X=TE5,Y=MASKEL)
      ENDIF     
C
C PREMIERS CALCULS POUR LES BILANS
C
      IF(BILMAS) THEN
C
        MASSES = 0.D0
        FLUSOR = 0.D0
        FLUENT = 0.D0
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE BILAN'
        CALL BILAN(MESH,H,T1,MASK,AT,0.D0,LT,NIT,LISTIN,
     *             MASSES,MSK,MASKEL,EQUA,TE5,OPTBAN,
     *             MESH%NPTFR,FLBOR,
     *             FLUX_BOUNDARIES,NUMLIQ%I,NFRLIQ)
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE BILAN'
C
        IF(NTRAC.GT.0) THEN
C
          IF(EQUA(1:15).NE.'SAINT-VENANT VF') THEN
            IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE BILANT'
            DO ITRAC=1,NTRAC
            MASSOU(ITRAC) = 0.D0
            CALL BILANT(H,T2,T3,DT,LT,NIT,LISTIN,
     *                  T%ADR(ITRAC)%P,
     *                  AGGLOT,MASSOU(ITRAC),MASTR0(ITRAC),
     *                  MASTR2(ITRAC),MASTEN(ITRAC),
     *                  MASTOU(ITRAC),MSK,MASKEL,MESH,FLBOR,
     *                  NUMLIQ%I,NFRLIQ,NPTFR,NAMETRAC(ITRAC),
     *                  FLBORTRA)
            ENDDO
            IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE BILANT'
C
          ELSE
            FLUTSOR = 0.D0
            FLUTENT = 0.D0
            IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE BILANT1'
            DO ITRAC=1,NTRAC
            CALL BILANT1(H,UCONV,VCONV,HPROP,T2,T3,T4,T5,T6,
     *                   DT,LT,NIT,ENTET,MASKTR,
     *                   T%ADR(1)%P,TN%ADR(1)%P,TETAT,
     *                   MASSOU(ITRAC),MSK,MASKEL,MESH,
     *                   FLUTSOR(ITRAC),FLUTENT(ITRAC),EQUA,LTT,ITRAC)
            ENDDO
            IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE BILANT1'
          ENDIF
C
        ENDIF
C
      ENDIF
C
C=======================================================================
C
      IF(NIT.EQ.0) THEN
        IF(LISTIN) THEN
          IF(LNG.EQ.1) WRITE(LU,9)
          IF(LNG.EQ.2) WRITE(LU,10)
        ENDIF
9      FORMAT(1X,'ARRET DANS TELEMAC, NOMBRE D''ITERATIONS DEMANDE NUL')
10     FORMAT(1X,'STOP IN TELEMAC, NUMBER OF TIME STEP ASKED EQUALS 0')
       STOP
      ENDIF
C
C=======================================================================
C
C
C     CAS D'UN COUPLAGE AVEC SISYPHE
C     ECRITURE DES CONDITIONS INITIALES DE U, V ET H
C
      IF(COUPLING.NE.' ') THEN
        IF(LNG.EQ.1) WRITE(LU,*) 'TELEMAC2D COUPLE AVEC : ',COUPLING
        IF(LNG.EQ.2) WRITE(LU,*) 'TELEMAC2D COUPLED WITH: ',COUPLING
      ENDIF
C
C     INITIALISATION DES VARIABLES CONSTANT FLOW DISCHARGE (CF. SISYPHE)
C     ------------------------------------------------------------------
C
      SISYPHE_CFD   = .FALSE.
      CONSTFLOW_SIS = .FALSE.
      NSIS_CFD      = 1
C
      IF(INCLUS(COUPLING,'SISYPHE')) THEN
C
         IF(INCLUS(COUPLING,'FILE-SISYPHE')) THEN
C
           WRITE (LU,*) 'TELEMAC-2D: FILE-COUPLING HAS NOW BEEN'
           WRITE (LU,*) '            SUPPRESSED'
           WRITE (LU,*) '            USE INTER-SISYPHE OR SISYPHE'
           WRITE (LU,*) '            INSTEAD OF FILE-SISYPHE'
           CALL PLANTE(1)
           STOP
C
         ELSEIF(INCLUS(COUPLING,'SISYPHE')) THEN
C
           IF(LNG.EQ.1) THEN
             WRITE (LU,*) 'TELEMAC-2D : COUPLAGE INTERNE AVEC SISYPHE'
           ENDIF
           IF(LNG.EQ.2) THEN
             WRITE (LU,*) 'TELEMAC-2D: INTERNAL COUPLING WITH SISYPHE'
           ENDIF
           CALL CONFIG_CODE(2) 
           IF(DEBUG.GT.0) WRITE(LU,*) 'PREMIER APPEL DE SISYPHE'
           CALL SISYPHE(0,LT,LEOPRD,LISPRD,NIT,U,V,H,H,ZF,CF,CF,
     *                  CONSTFLOW_SIS,NSIS_CFD,SISYPHE_CFD,CODE1,PERCOU,
     *                  U,V,AT,VISC,DT,CHARR,SUSP,
C                                      CHARR,SUSP : RETURNED BY SISYPHE
C                                                   BUT THEN GIVEN TO IT
     *                  FLBOR,SOLSYS,DM1,USIS,VSIS,ZCONV)
           IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE SISYPHE'
           CALL CONFIG_CODE(1)
C          AVOIDING TWO OUTPUTS WHEN SISYPHE IS CALLED TWICE
           IF(SUSP.AND.CHARR.AND.PERCOU.NE.1) THEN
             LEOPRD_CHARR=NIT+PERCOU
           ELSE
             LEOPRD_CHARR=LEOPRD
           ENDIF   
C
         ENDIF

       ENDIF
C
C=======================================================================
C INITIALISE INFILTRATION STRUCTURES FOR THE COUPLING WITH ESTEL3D
C
      CALL INFILTRATION_INIT(NPOIN,(CODE(1:7).EQ.'ESTEL3D'))
C
C     SAVE THE DEPTH CALCULATED BY TELEMAC2D FOR ESTEL3D
C
      IF(CODE(1:7).EQ.'ESTEL3D') CALL DEPTH_FILL(H%R)
C
C=======================================================================
C
C : 3                    /* BOUCLE EN TEMPS */
C
C=======================================================================
C
C     STORAGE OF DT FOR CASE WITH VARIABLE TIME-STEP
C
      DTCAS = DT
C
C     CALLED BY ANOTHER PROGRAM, ONLY INITIALISATION REQUIRED
      IF(PASS.EQ.0) THEN
        IF(LNG.EQ.1) WRITE(LU,*) 'FIN D''INITIALISATION DE TELEMAC2D'
        IF(LNG.EQ.2) WRITE(LU,*) 'TELEMAC2D INITIALISED'
        RETURN
      ENDIF
C
700   CONTINUE
C
      IF(PASS.EQ.1) THEN
        IF(CODE(1:7).EQ.'ESTEL3D') THEN
          AT=ATDEP
          NIT=NITER
! --- JP RENAUD START ---
C         USE THE TIME STEP SPECIFIED BY ESTEL-3D
          IF(PRESENT(DTDEP)) THEN
            DT = DTDEP
            DTCAS = DTDEP
          ! TODO: CHECK WHAT HAPPENS WITH ADAPTIVE TIME STEP
          ENDIF
! --- JP RENAUD END ---
        ELSE
          CALL PLANTE(1)
          STOP 'UNKNOWN CALLING PROGRAM'
        ENDIF
      ENDIF
C
      LT = LT + 1
C
      IF(DTVARI.AND.EQUA(1:15).NE.'SAINT-VENANT VF') THEN
C       NOMBRE DE COURANT SELON SCHEMA PSI EN P1
        CALL CFLPSI(T1,U,V,DT,IELM,MESH,MSK,MASKEL)
        CALL MAXI(CFLMAX,IMAX,T1%R,NPOIN)
C       LIMITATION A DES VARIATIONS DANS L'INTERVALLE ( 1/2 , 2 )
        DT = DT * MAX(MIN(CFLWTD/MAX(CFLMAX,1.D-6),2.D0),0.5D0)
C       LIMITATION AU DT DU FICHIER CAS
        DT=MIN(DT,DTCAS)
        IF(NCSIZE.GT.1) DT=P_DMIN(DT)
        IF(ENTET) THEN
          IF (LNG.EQ.1) WRITE(LU,78) CFLMAX,DT
          IF (LNG.EQ.2) WRITE(LU,79) CFLMAX,DT
78        FORMAT(1X,'    NOMBRE DE COURANT MAXIMUM :',G16.7,/,1X,
     *              '    PAS DE TEMPS              :',G16.7)
79        FORMAT(1X,'    MAXIMUM COURANT NUMBER: ',G16.7,/,1X,
     *              '    TIME-STEP                 :',G16.7)
        ENDIF
      ENDIF
C
C=======================================================================
C
C     COUPLAGE AVEC SISYPHE
C
      IF(INCLUS(COUPLING,'SISYPHE')) THEN
C
        CALL CONFIG_CODE(2)
C       HN NOT DEFINED HERE AT FIRST ITERATION
        IF(LT.EQ.1) CALL OS('X=Y     ',X=HN,Y=H)
C
        SUSP1=SUSP.AND.PERCOU.EQ.1
        IF(SUSP1.OR.(CHARR.AND.(PERCOU*((LT-1)/PERCOU).EQ.LT-1))) THEN
C     
          IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE SISYPHE, CHARRIAGE'        
          CALL SISYPHE(1,LT,LEOPRD_CHARR,LISPRD,NIT,U,V,H,HN,ZF,
     *                 CF,CF,CONSTFLOW_SIS,NSIS_CFD,SISYPHE_CFD,CODE1,
     *                 PERCOU,U,V,AT,VISC,DT*PERCOU,CHARR,SUSP1,
     *                 FLBOR,SOLSYS,DM1,USIS,VSIS,ZCONV)
          IF(DEBUG.GT.0) WRITE(LU,*) 'FIN APPEL SISYPHE, CHARRIAGE'
C             
        ENDIF
C
        IF(SUSP.AND.PERCOU.NE.1) THEN
C
          IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE SISYPHE, SUSPENSION'
          CALL SISYPHE(1,LT,LEOPRD,LISPRD,NIT,U,V,H,HN,ZF,
     *                 CF,CF,CONSTFLOW_SIS,NSIS_CFD,SISYPHE_CFD,CODE1,
     *                 1,U,V,AT,VISC,
     *                 DT,NON,SUSP,
     *                 FLBOR,SOLSYS,DM1,USIS,VSIS,ZCONV)    
          IF(DEBUG.GT.0) WRITE(LU,*) 'FIN APPEL DE SISYPHE, SUSPENSION'
C
        ENDIF
C 
        CALL CONFIG_CODE(1)
C  
      ENDIF
C
C=======================================================================
C
      IF(ADJO) THEN
        AT = AT - DT
      ELSE
C       DT IS NOT YET KNOWN IN FINITE VOLUMES
        IF(EQUA(1:15).NE.'SAINT-VENANT VF') AT = AT + DT
      ENDIF
C
      IF(DTVARI) THEN
        IF(AT.GT.DUREE+AT0) THEN
C         DERNIER PAS DE TEMPS
          NIT = LT
        ELSE
C         VALEUR BIDON SUPERIEURE A LT
          NIT = LT + 10
        ENDIF
      ENDIF
C
      IF((LISPRD*(LT/LISPRD).EQ.LT.AND.LT.GE.PTINIL).OR.LT.EQ.NIT) THEN
        ENTET=LISTIN
      ELSE
        ENTET=.FALSE.
      ENDIF
      
! --- JP RENAUD START ---
C CONSTRAIN TELEMAC-2D OUTPUT IN THE LISTING
      IF (PRESENT(DOPRINT)) ENTET = ENTET .AND. DOPRINT
! --- JP RENAUD END ---

      IF(ENTET) CALL ENTETE(2,AT,LT)
C
C=======================================================================
C
C SAUVEGARDE DE UN, VN, HN, TN, AKN ET EPN (ILS SONT DANS LE BLOC FN)
C
! --- JP RENAUD START ---
C THIS IS NOT DONE WHEN ITERATING FOR THE COUPLING WITH ESTEL-3D
       IF(CODE(1:7).EQ.'ESTEL3D'.AND.PRESENT(NEWTIME)) THEN
        IF(NEWTIME) CALL OS('X=Y     ',X=FN,Y=F)
       ELSE
         CALL OS('X=Y     ',X=FN,Y=F)
       ENDIF
!      CALL OS( 'X=Y     ' , FN , F , F , C )
! --- JP RENAUD END ---     
C
C=======================================================================
C
C NEW COUPLING WITH SISYPHE FOR CONSTANT FLOW DISCHARGE
C
      IF(SISYPHE_CFD.AND.CONSTFLOW_SIS) GOTO 999
C
      DO 888 ISIS_CFD=1,NSIS_CFD
C
C=======================================================================
C
C  MASQUAGE DES ELEMENTS DECOUVRANTS
C
      IF(MSK) CALL OS( 'X=C     ' , MASKEL , S , S , 1.D0 )
      IF (OPTBAN.EQ.2) THEN
        CALL MASKBD(MASKEL%R,ZFE%R,ZF%R,H%R,
     *              HMIN,MESH%IKLE%I,MESH%IFABOR%I,IT1%I,NELEM,NPOIN)
      ENDIF
C
C  MASQUAGE ARTIFICIEL DONNE PAR L'UTILISATEUR
C
      IF(MSKUSE) THEN
      CALL MASKOB (MASKEL%R,MESH%X%R,MESH%Y%R,
     *             MESH%IKLE%I,NELEM,NELMAX,NPOIN,AT,LT)
      ENDIF
C
C CONSTRUCTION DU MASQUE DES POINTS A PARTIR DU MASQUE DES ELEMENTS
C ET CHANGEMENT DE IFAMAS (IFABOR AVEC MASQUAGE)
C
      IF(MSK) THEN
        CALL MASKTO(MASKEL%R,MASKPT,IFAMAS%I,MESH%IKLE%I,
     *              MESH%IFABOR%I,MESH%ELTSEG%I,MESH%NSEG,
     *              NELEM,NPOIN,IELMT,MESH)
        IF(IELMX.NE.IELM1) CALL CHGDIS(MASKPT,IELM1,IELMX,MESH)
      ENDIF
C
C-----------------------------------------------------------------------
C  CALCUL DE L'INTEGRALE DES BASES
C-----------------------------------------------------------------------
C
!     IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE MASBAS2D'
!     IF(MSK) THEN
!       CALL MASBAS2D(VOLU2D,V2DPAR,UNSV2D,IELM1,MESH,MSK,MASKEL,T2,T2)
!     ENDIF
!     IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE MASBAS2D'
C
C-----------------------------------------------------------------------
C
C UPDATING POROSITY : NEW VALUE IN TE5
C                     OLD - NEW IN TE4
C
      IF(OPTBAN.EQ.3) THEN
C
         CALL OS('X=Y     ',TE4,TE5,TE5,0.D0)
         CALL POROS(TE5,ZF,HN,MESH)
         IF(MSK) CALL OS('X=XY    ',X=TE5,Y=MASKEL)                      
C        TEST OF UNDER-RELAXATION
         RELAX = 0.05D0
         CALL OS('X=CX    ',X=TE5,C=RELAX)
         CALL OS('X=X+CY  ',X=TE5,Y=TE4,C=1.D0-RELAX)
C        TE4 = OLD POROS - NEW POROS
         CALL OS('X=X-Y   ',X=TE4,Y=TE5)
C
      ENDIF
C
C=======================================================================
C
C NOUVEAUX CHAMPS CONVECTEUR ET PROPAGATEUR
C REMARQUER QU'A CE NIVEAU, U = UN  V = VN ET H = HN
C
      IF(CONV) THEN
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE CHPCON'
        CALL CHPCON(UCONV,VCONV,U,V,UN,VN,TETAU)
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE CHPCON'
      ENDIF
C
C     CALCUL DU NOUVEAU PROPAGATEUR
C
      IF(PROPA) THEN
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE HPROPA'
        CALL HPROPA(HPROP ,HN,H,PROLIN,HAULIN,TETAHC,NSOUSI)
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE HPROPA'
      ENDIF
C
C=======================================================================
C
C PREPARATION DES CONDITIONS AUX LIMITES POUR LES SEUILS.
C
      IF(NWEIRS.GT.0) THEN
C
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE CLSING'
        CALL CLSING(NWEIRS,NPSING,NPSMAX,NUMDIG%I,
     *              MESH%X%R,MESH%Y%R,ZF%R,CHESTR%R,NKFROT%I,
     *              KARMAN,ZDIG%R,PHIDIG%R,MESH%NBOR%I,
     *              H%R,T,NTRAC,IOPTAN,T1%R,UBOR%R,VBOR%R,TBOR,
     *              LIHBOR%I,LIUBOR%I,LIVBOR%I,LITBOR,GRAV)
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE CLSING'
C
      ENDIF
C
C ON CONSIDERE QUE LES TYPES DE CONDITIONS AUX LIMITES NE CHANGENT PAS
C PENDANT LES SOUS-ITERATIONS. SINON IL FAUT DEPLACER LES APPELS DE
C KEPSIN,DIFFIN,PROPIN
C
C TYPES DE CONDITIONS POUR LE MODELE K-EPSILON
C
      IF(ITURB.EQ.3) CALL KEPSIN(LIMKEP%I,LIUBOR%I,NPTFR,
     *                           KENT,KENTU,KSORT,KADH,KLOG,
     *                           KINC,KNEU,KDIR)
C
C TYPES DE CONDITIONS POUR LA DIFFUSION DU TRACEUR :
C
      IF(NTRAC.GT.0) THEN
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE DIFFIN'
        DO ITRAC=1,NTRAC 
        CALL DIFFIN(MASKTR,LIMTRA%I,LITBOR%ADR(ITRAC)%P%I,
     *              IT1%I,U%R,V%R,MESH%XNEBOR%R,MESH%YNEBOR%R,
     *              MESH%NBOR%I,MESH%KP1BOR%I,NPTFR,
     *              KENT,KSORT,KLOG,KINC,KNEU,KDIR,KDDL,
     *              ICONVF(3),MESH%NELBOR%I,NPOIN,NELMAX,MSK,MASKEL%R,
     *              NFRLIQ,THOMFR,DEBLIQ,FINLIQ,FRTYPE,
     *              TN%ADR(ITRAC)%P,TBOR%ADR(ITRAC)%P,MESH)
        ENDDO
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE DIFFIN'
      ENDIF
C
C TYPES DE CONDITIONS POUR LA PROPAGATION:
C NECESSAIRE POUR THOMFR ?? (SINON REFAIT APRES BORD !)
C
C
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE PROPIN'
      CALL PROPIN_TELEMAC2D 
     *            (LIMPRO%I,LIMPRO%DIM1,MASK,LIUBOR%I,LIVBOR%I,
     *             LIHBOR%I,MESH%KP1BOR%I,MESH%NBOR%I,NPTFR,
     *             KENT,KENTU,KSORT,KADH,KLOG,KINC,
     *             KNEU,KDIR,KDDL,KOND,CLH%I,CLU%I,CLV%I,
     *             U%ELM,U%R,V%R,GRAV,H%R,LT,NPOIN,
     *             MESH%NELBOR%I,NELMAX,MSK,MASKEL%R,
     *             NFRLIQ,THOMFR,DEBLIQ,FINLIQ,FRTYPE,
     *             MESH%XNEBOR%R,MESH%YNEBOR%R,.FALSE.,MESH)
C    *             MESH%XNEBOR%R,MESH%YNEBOR%R, ENTET ,MESH)
C       WARNINGS WILL BE GIVEN AT THE SECOND CALL AFTER BORD
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE PROPIN'
C
C=======================================================================
C                 CALCUL DES COEFFICIENTS DE FROTTEMENT
C                         VARIABLES EN TEMPS
C=======================================================================
C CORSTR NE FAIT RIEN SAUF SI IL EST ECRIT PAR L'UTILISATEUR.
C
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE CORSTR'
      CALL CORSTR
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE CORSTR'
C
      IF(IELMU.EQ.12.OR.IELMU.EQ.13) CALL CHGDIS(CHESTR,11,IELMU,MESH)
C
      CALL FRICTION_CHOICE(1,KARMAN)
C
C=======================================================================
C                 CALCUL DES COEFFICIENTS DE VISCOSITE
C=======================================================================
C
C  CALCUL DE LA VISCOSITE DYNAMIQUE VISC
C
      IF(ITURB.EQ.1) THEN
C
        CALL OS( 'X=C     ' , VISC , VISC , VISC , PROPNU )
C
      ELSEIF(ITURB.EQ.2) THEN
C
        CALL DISPER( VISC , U%R , V%R , H%R , CF%R , ELDER , PROPNU )
C
      ELSEIF(ITURB.EQ.3) THEN
C
        CALL VISTUR(VISC,AK,EP,NPOIN,CMU,PROPNU)
C
      ELSEIF(ITURB.EQ.4) THEN
C
        CALL SMAGOR(VISC,CF,U,V,MESH,T1,T2,T3,T4,MSK,MASKEL,PROPNU)
C
      ELSE
C
        IF(LISTIN) THEN
          IF(LNG.EQ.1) WRITE(LU,15) ITURB
          IF(LNG.EQ.2) WRITE(LU,16) ITURB
        ENDIF
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C  COEFFICIENT DE DIFFUSION THERMIQUE (PRANDTL = 1 POUR L'INSTANT)
C  ET LE MEME POUR TOUS LES TRACEURS
C
      IF(NTRAC.GT.0.AND.DIFT) THEN
        DO ITRAC=1,NTRAC
          CALL OS( 'X=Y     ' , X=VISCT%ADR(ITRAC)%P , Y=VISC )
          CALL OS( 'X=X+C   ' , X=VISCT%ADR(ITRAC)%P , C=DIFNU-PROPNU )
        ENDDO
      ENDIF
C
C  CORRECTION EVENTUELLE DES COEFFICIENTS DE VISCOSITE.
C
      CALL CORVIS
!
!=======================================================================
!  SOURCES : COMPUTATION OF INPUTS WHEN VARYING IN TIME
!            IF NO VARIATION IN TIME DSCE2=DSCE AND TSCE2=TSCE
!=======================================================================
!
      IF(NREJET.GT.0) THEN
        DO I=1,NREJET
          DSCE2(I)=DEBSCE(AT,I,DSCE)
        ENDDO
        IF(NTRAC.GT.0) THEN
          DO I=1,NREJET 
            DO ITRAC=1,NTRAC                       
              TSCE2(I,ITRAC)=TRSCE(AT,I,ITRAC)
            ENDDO
          ENDDO
        ENDIF
      ENDIF
C
C=======================================================================
C CONDITIONS AUX LIMITES
C=======================================================================
C
      IF(THOMFR) THEN
C
      CALL CPSTVC(H,T9)
      CALL PREBOR(HBOR%R,UBOR%R,VBOR%R,TBOR,U%R,V%R,H%R,
     *            T9%R,T,MESH%NBOR%I,MESH%KP1BOR%I,
     *            NPOIN,NPTFR,NTRAC,DEBLIQ,FINLIQ,NFRLIQ)
C
      ENDIF
C
C APPEL DU SOUS-PROGRAMME UTILISATEUR POUR LES CONDITIONS AUX LIMITES.
C
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE BORD'
      CALL BORD(HBOR%R,UBOR%R,VBOR%R,TBOR,
     *          U,V,H,ZF%R,MESH%NBOR%I,W1,T8,
     *          LIHBOR%I,LIUBOR%I,LITBOR,
     *          MESH%XNEBOR%R,MESH%YNEBOR%R,NPOIN,NPTFR,
     *          NPTFR2,AT,
     *          NDEBIT,NCOTE,NVITES,NTRAC,NTRACE,NFRLIQ,NUMLIQ%I,
     *          KENT,KENTU,PROVEL,MASK,MESH,EQUA,T2D_FILES(T2DIMP)%NAME)
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE BORD'
C
C CALCUL DE LIMPRO, CLU,CLV, CLH ET MASK
C    
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE PROPIN'
      CALL PROPIN_TELEMAC2D 
     *            (LIMPRO%I,LIMPRO%DIM1,MASK,LIUBOR%I,LIVBOR%I,
     *             LIHBOR%I,MESH%KP1BOR%I,MESH%NBOR%I,NPTFR,
     *             KENT,KENTU,KSORT,KADH,KLOG,KINC,
     *             KNEU,KDIR,KDDL,KOND,CLH%I,CLU%I,CLV%I,
     *             U%ELM,U%R,V%R,GRAV,H%R,LT,NPOIN,
     *             MESH%NELBOR%I,NELMAX,MSK,MASKEL%R,
     *             NFRLIQ,THOMFR,DEBLIQ,FINLIQ,FRTYPE,
     *             MESH%XNEBOR%R,MESH%YNEBOR%R,ENTET,MESH)
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE PROPIN' 
C
C CONDITIONS AUX LIMITES EN K-EPSILON: KBOR,EBOR ET AUBOR
C
      IF(ITURB.EQ.3) THEN
         CALL KEPSCL(KBOR%R,EBOR%R,AUBOR%R,CF%R,CFBOR%R,
     *          MESH%DISBOR%R,
     *          UN%R,VN%R,HN%R,LIMKEP%I,LIUBOR%I,LIMPRO%I,
     *          MESH%NBOR%I,NPTFR,KARMAN,CMU,C2,ESTAR,SCHMIT,LISRUG,
     *          PROPNU,KMIN,EMIN,KNEU,KDIR,KENT,KENTU,KADH,KLOG)
      ENDIF  
C
C APPEL DU SYSTEME DE RESOLUTION AUX LIMITES PAR CARACTERISTIQUES
C (METHODE DE THOMPSON)
C
      IF(THOMFR) THEN
C
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE THOMPS'
      CALL THOMPS(HBOR%R,UBOR%R,VBOR%R,TBOR,U,V,T9,
     *            T,ZF,MESH%X%R,MESH%Y%R,MESH%NBOR%I,
     *            FRTYPE,T1,T2,T3,T4,T6,FU,FV,
     *            LIHBOR%I,LIUBOR%I,LIVBOR%I,LITBOR,IT1%I,
     *            T8,W1,IT2%I,
     *            CV2%R,CV3%R,TE1%R,TE2%R,HTILD%R,UTILD%R,VTILD%R,
     *            TTILD,TE3%R,
     *            MESH%SURDET%R,MESH%IKLE%I,CF,SMH,
     *            MESH%IFABOR%I,MESH%NULONE%I,NELEM,MESH,
     *            MESH%KP1BOR%I,MESH%XNEBOR%R,MESH%YNEBOR%R,
     *            NPOIN,NPTFR,
     *            LT,NIT,AT,DT,GRAV,DEBLIQ,FINLIQ,NTRAC,
     *            NFRLIQ,KSORT,MESH%LV,
     *            MSK,MASKEL,MASKPT,MESH%NELBOR%I,
     *            NELMAX,IELM,NORD,FAIR,WINDX,WINDY,
     *            VENT,HWIND,
     *            CORIOL,FCOR,SPHERI,OPTPRO,MAREE,MARDAT,MARTIM,
     *            PHI0,OPTSOU,ISCE,DSCE2,USCE,VSCE,T5%R,COUROU,NPTH,
     *            VARCL,NVARCL,VARCLA,NUMLIQ%I,BM1%X%R,UNSV2D,HFROT)
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE THOMPS'
C
      ENDIF
C
C     ON CONTROLE HBOR CAR L'UTILISATEUR PEUT RECRIRE BORD ET SE TROMPER
      CALL CLIP(HBOR,0.D0,.TRUE.,1.D6,.FALSE.,0)
C
C=======================================================================
C
C BOUCLE DES SOUS-ITERATIONS OU L'ON ACTUALISE CONVECTEUR ET PROPAGATEUR
C
C=======================================================================
C
      DO 701 ISOUSI = 1 , NSOUSI
      IF(DEBUG.GT.0) WRITE(LU,*) 'BOUCLE 701 ISOUSI=',ISOUSI
C
C=======================================================================
C
C : 4                     CONVECTION
C
C=======================================================================
C
      IF(CONV.AND.FTILD%N.GT.0) THEN
C
        IF(ENTET) CALL ENTETE(3,AT,LT)
C
        IF(SPHERI) THEN
          CALL OS('X=Y/Z   ',UCONV,UCONV,MESH%COSLAT,C)
          CALL OS('X=Y/Z   ',VCONV,VCONV,MESH%COSLAT,C)
        ENDIF
C
C       APPEL DE CHARAC
C
        CALL CHARAC( FNCAR , FTILD  , FTILD%N  , UCONV , VCONV,S,S,
     *               DT    , IFAMAS , IELM     , NPOIN , 1 , 1,
     *               MSK   , MASKEL , BM1%X    , BM1%D , TB   ,
     *               IT1%I , IT2%I  , IT3%I    , IT4%I , MESH ,
     *               MESH%NELEM,MESH%NELMAX,MESH%IKLE,MESH%SURDET)
C
        IF(SPHERI) THEN
          CALL OS('X=XY    ',UCONV,MESH%COSLAT,S,C)
          CALL OS('X=XY    ',VCONV,MESH%COSLAT,S,C)
        ENDIF
C
        IF(IELM1.NE.IELMH.AND.CONVV(2)) THEN
          CALL CHGDIS(HTILD,IELM1,IELMH,MESH)
        ENDIF
        IF(IELM1.NE.IELMU.AND.CONVV(1)) THEN
C         POINTS QUASI-BULLE OBTENUS PAR INTERPOLATION
          IF(IELMU.EQ.12) THEN
            CALL CHGDIS(UTILD,IELM1,IELMU,MESH)
            CALL CHGDIS(VTILD,IELM1,IELMU,MESH)
          ENDIF
        ENDIF
        IF(NTRAC.GT.0.AND.IELM1.NE.IELMT.AND.CONVV(3)) THEN
          DO ITRAC=1,NTRAC
            CALL CHGDIS(TTILD%ADR(ITRAC)%P,IELM1,IELMT,MESH)
          ENDDO
        ENDIF
C
      ENDIF
C
C GESTION DES TABLEAUX .
C
      CALL GESTIO(UN   ,VN   ,HN   ,TN   ,AKN   ,EPN   ,
     *            UTILD,VTILD,HTILD,TTILD,AKTILD,EPTILD,
     *            NTRAC.GT.0,PROPA,CONVV,ITURB,3)
C
C=======================================================================
C                  FIN DE L'ETAPE DE CONVECTION
C=======================================================================
C=======================================================================
C
C : 6                DIFFUSION - PROPAGATION
C
C=======================================================================
C
      IF(PROPA) THEN
      IF(ENTET) CALL ENTETE(6,AT,LT)
C     L'INFORMATION SUR LA METHODE DE RESOLUTION N'EST DONNEE QU'EN CAS
C     DE SORTIE LISTING
      INFOGS=.FALSE.
      IF(INFOGR.AND.ENTET) INFOGS=.TRUE.
C
C  CONDITIONS METEO.
C
      IF(VENT.OR.ATMOS) THEN
        CALL METEO(PATMOS%R,WINDX%R,WINDY%R,
     *             FUAIR,FVAIR,MESH%X%R,MESH%Y%R,AT,LT,NPOIN,VENT,ATMOS,
     *             H%R,T1%R,GRAV,ROEAU,NORD,PRIVE)
      ENDIF
C
C  CALCUL DE LA MASSE VOLUMIQUE QUAND ELLE EST VARIABLE
C
      IF(ROVAR) THEN
C       BEWARE, SALINITY MUST BE HERE THE FIRST TRACER
        CALL VALRO(RO,T,ROEAU)
      ENDIF
C
C  TERMES SOURCES DUS AUX BUSES ET SIPHONS.
C
      IF(NSIPH.GT.0) THEN
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE SIPHON'
      CALL SIPHON(RELAXS,NSIPH,ENTSIP,SORSIP,GRAV,
     *            H%R,ZF%R,ISCE,DSCE,SECSCE,ALTSCE,CSSCE,CESCE,
     *            DELSCE,ANGSCE,LSCE,
     *            NTRAC,T,TSCE,USCE,VSCE,U%R,V%R,ENTET,MAXSCE)
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE SIPHON'
      ENDIF
C
C  TERMES SOURCES DE L'ETAPE DE PROPAGATION .
C
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE PROSOU'
      CALL PROSOU(FU,FV,SMH,UN,VN,HN,GRAV,NORD,
     *            FAIR,WINDX,WINDY,VENT,HWIND,
     *            CORIOL,FCOR,SPHERI,YASMH,
     *            MESH%COSLAT,MESH%SINLAT,AT,LT,
     *            NREJET,NREJEU,DSCE2,ISCE,T1,MESH,MSK,MASKEL,
     *            MAREE,MARDAT,MARTIM,PHI0,OPTSOU,COUROU,NPTH,
     *            VARCL,NVARCL,VARCLA,UNSV2D)
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE PROSOU'
C
C  ETAPE DE PROPAGATION.
C
      IF(EQUA(1:15).EQ.'SAINT-VENANT EF'.OR.
     *   EQUA(1:10).EQ.'BOUSSINESQ') THEN
C
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE PROPAG'
      CALL PROPAG
     *(U,V,H,UCONV,VCONV,CONVV,H0,C0,COTOND,PATMOS,ATMOS,
     * HPROP,UN,VN,HN,UTILD,VTILD,HTILD,DH,DU,DV,DHN,VISC,VISC_S,
     * FU,FV,
     * SMH,MESH,ZF,AM1,AM2,AM3,BM1,BM2,CM1,CM2,TM1,A23,A32,MBOR,
     * CV1,CV2,CV3,W1,UBOR,VBOR,AUBOR,HBOR,DIRBOR,
     * TE1,TE2,TE3,TE4,TE5,T1,T2,T3,T4,T5,T6,T7,T8,
     * LIMPRO,MASK,GRAV,ROEAU,CF,DIFVIT,IORDRH,IORDRU,LT,AT,DT,
     * TETAC,TETAHC,TETAU,TETAD,
     * AGGLOC,AGGLOU,KDIR,INFOGS,KFROT,ICONVF,
     * PRIVE,ISOUSI,BILMAS,MASSES,YASMH,OPTBAN,CORCON,
     * OPTSUP,MSK,MASKEL,MASKPT,RO,ROVAR,
     * MAT,RHS,UNK,TB,S,TB,PRECCU,SOLSYS,CFLMAX,OPDVIT,
C                       TB HERE TO REPLACE BD SUPPRESSED, NOT USED
     * OPTSOU,NFRLIQ,SLVPRO,EQUA,VERTIC,ADJO,ZFLATS,TETAZCOMP,
     * UDEL,VDEL,DM1,ZCONV,COUPLING,FLBOR,BM1S,BM2S,CV1S,
     * VOLU2D,V2DPAR,UNSV2D,NUMDIG%I,NWEIRS,NPSING,HFROT)
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE PROPAG'
C
      IF(ADJO) THEN
C
       IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE PROPAG_ADJ'
       CALL PROPAG_ADJ
     *(UCONV,VCONV,CONVV,H0,C0,COTOND,PATMOS,ATMOS,
     * HPROP,UN,VN,HN,UTILD,VTILD,HTILD,DH,DU,DV,DHN,VISC,VISC_S,
     * FU,FV,SMH,MESH,ZF,AM1,AM2,AM3,BM1,BM2,CM1,CM2,TM1,A23,A32,
     * MBOR,CV1,CV2,CV3,W1,UBOR,VBOR,AUBOR,HBOR,DIRBOR,
     * TE1,TE2,TE3,TE4,TE5,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,
     * LIMPRO,MASK,GRAV,ROEAU,CF,DIFVIT,IORDRH,IORDRU,LT,AT,DT,
     * TETAC,TETAHC,TETAU,TETAD,
     * AGGLOC,AGGLOU,KDIR,INFOGS,KFROT,ICONVF,
     * PRIVE,ISOUSI,BILMAS,MASSES,YASMH,OPTBAN,CORCON,
     * OPTSUP,MSK,MASKEL,MASKPT,RO,ROVAR,
     * MAT,RHS,UNK,TB,S,TB,PRECCU,SOLSYS,CFLMAX,OPDVIT,
     * OPTSOU,NFRLIQ,SLVPRO,EQUA,VERTIC,
     * ADJO,UD,VD,HD,U,V,H,UU,VV,HH,UIT1,VIT1,HIT1,PP,QQ,RR,
     * TAM1,TAM2,TAM3,TBM1,TBM2,TCM1,TCM2,MATADJ,UNKADJ,
     * ALPHA1,ALPHA2,ALPHA3,ADJDIR,ESTIME,OPTCOST,NIT,NVARRES,
     * VARSOR,T2D_FILES(T2DRES)%LU,T2D_FILES(T2DREF)%LU,
     * ALIRE,TROUVE,MAXVAR,VARCL,VARCLA,TEXTE,
     * TEXREF,TEXRES,W,OUTINI,CHESTR,KARMAN,NDEF,ITURB,LISRUG,
     * LINDNER,SB,DP,SP,CHBORD,CFBOR,HFROT,UNSV2D)
       IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE PROPAG_ADJ'
C
      ENDIF 
C
      ELSEIF(EQUA(1:15).EQ.'SAINT-VENANT VF') THEN
C
C      VOLFIN MAY CHANGE DT
C
C      CM1%D%R : HT  
C
       IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE VOLFIN'
       CALL VOLFIN(W1%R,AT,DT,LT,NIT,NELEM,NPTFR,
     *      TB,ZF%R,CHESTR%R,NPOIN,HN%R,H%R,U%R,V%R,QU%R,QV%R,
     *      GRAV,ENTET,S,MSK,MASKEL,MESH,LIMPRO%I,
     *      MESH%NBOR%I,KDIR,KNEU,KDDL,HBOR%R,UBOR%R,VBOR%R,
     *      MASSES,FLUENT,FLUSOR,CFLWTD,DTVARI,KFROT,
     *      NREJET,ISCE,TSCE2,MAXSCE,MAXTRA,YASMH,SMH%R,
     *      NTRAC,T%ADR(1)%P%DIM1,T,HT,TN,
     *      LIMTRA%DIM1,LIMTRA%I,
     *      TBOR,MASSOU,FLUTENT,FLUTSOR,MESH%DTHAUT%R,
     *      MESH%DPX%R,MESH%DPY%R,CM1%X%R,CM2%X%R,
     *      MESH%CMI%R,MESH%JMI%I,TE1%R,TE2%R,
     *      DIFVIT,ITURB,PROPNU,DIFT,DIFNU,
     *      BM1%X%R,BM2%X%R,OPTVF,
     *      HSTOK%R,HCSTOK%R,LOGFR%I,DSZ%R,FLUXT,FLUHBOR,
     *      DTN,FLUSORTN,FLUENTN,
     *      LTT,FLUXTEMP,FLUHBTEMP,HC%R,SMTR,MESH%AIRST%R,
     *      TMAX,DTT) 
       IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE VOLFIN'
C
       AT = AT + DT
       IF(AT.GE.TMAX) THEN 
         NIT = LT
         IF(LISTIN) CALL ENTETE(1,AT,LT)
       ENDIF
C
      ELSE
C
        IF(LNG.EQ.1) WRITE(LU,*) 'EQUATIONS INCONNUES : ',EQUA
        IF(LNG.EQ.2) WRITE(LU,*) 'UNKNOWN EQUATIONS: ',EQUA
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C SI IL N'Y A PAS D'ETAPE DE PROPAGATION :
C
      ELSE
C
C GESTION DES TABLEAUX .
C
         CALL GESTIO(U    ,V    ,H    ,T,AK  ,EP ,
     *               UTILD,VTILD,HTILD,T,AK  ,EP ,
     *               NTRAC.GT.0,PROPA,CONVV,ITURB ,6)
C
C        SMH UTILISE PAR LE TRACEUR
C        SI ON VEUT SIMULER SUBIEF EN ENLEVANT PROPAGATION
C        ET CONVECTION, PROSOU N'EST PAS APPELE ET LES SOURCES
C        PONCTUELLES NE SONT PAS PRISES EN COMPTE.
C        EN TOUTE RIGUEUR IL FAUDRAIT METTRE ICI CALL PROSOU.
         IF(NTRAC.GT.0) CALL OS('X=0     ',X=SMH)
C
      ENDIF
C
C     TREATMENT OF NEGATIVE DEPTHS
C
      CALL CORRECTION_DEPTH_2D(MESH%GLOSEG%I,MESH%GLOSEG%DIM1,
     *                         YAFLODEL,YASMH)
C
C=======================================================================
C                        FIN DE LA PROPAGATION
C=======================================================================
C
C  CALCUL DES NOUVEAUX CHAMPS CONVECTEURS SI IL RESTE DES
C  SOUS-ITERATIONS.
C
C  LE TEST SUR ISOUSI N'EST FAIT QUE POUR HPROP ET PAS POUR UCONV
C  POUR DES RAISONS DE CONSERVATION DE LA MASSE DE TRACEUR
C  (IL FAUT GARDER LE MEME HPROP POUR LE TRACEUR QUE POUR H ET U)
C
      IF(ISOUSI.NE.NSOUSI) THEN
C      CALCUL DU NOUVEAU PROPAGATEUR SI IL Y A PROPAGATION
       IF(PROPA) THEN
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE HPROPA'
        CALL HPROPA(HPROP ,HN,H,PROLIN,HAULIN,TETAHC,NSOUSI)
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE HPROPA'
       ENDIF
      ENDIF
C
C     CALCUL DU NOUVEAU CHAMP CONVECTEUR (S'IL Y A CONVECTION)
      IF(CONV) THEN
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE CHPCON'
        CALL CHPCON(UCONV,VCONV,U,V,UN,VN,TETAU)
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE CHPCON'
      ENDIF
C
C=======================================================================
C FIN DE LA BOUCLE DES SOUS-ITERATIONS
C
701   CONTINUE
C
C=======================================================================
C
C : 5                   DIFFUSION DU TRACEUR
C
C=======================================================================
C
      IF(NTRAC.GT.0.AND.EQUA(1:15).NE.'SAINT-VENANT VF') THEN
C
      IF(ENTET) CALL ENTETE(5,AT,LT)
C
      DO ITRAC=1,NTRAC
C
C  CONDITIONS AUX LIMITES POUR L'ETAPE DE DIFFUSION DU TRACEUR.
C
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE DIFFCL POUR ITRAC=',ITRAC
      CALL DIFFCL(LITBOR%ADR(ITRAC)%P%I,
     *            TTILD%ADR(ITRAC)%P%R,TBOR%ADR(ITRAC)%P%R,
     *            MESH%NBOR%I,ICONVF(3),NPOIN,NPTFR)
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE DIFFCL'
C
      ENDDO
C
C  TERMES SOURCES DE L'ETAPE DE DIFFUSION-TERMES SOURCES DU TRACEUR
C
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE DIFSOU'
      CALL DIFSOU(TEXP,TIMP,YASMI,TSCEXP,HPROP,TN,TETAT,NREJTR,
     *            ISCE,DSCE2,TSCE2,MAXSCE,MAXTRA,AT,DT,MASSOU,NTRAC,
     *            MESH%FAC%R)
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE DIFSOU'
C
      DO ITRAC=1,NTRAC
C
C  APPEL DU DIFFUSEUR STANDARD. (CV1 EST LE SECOND MEMBRE)
C
      INFOGT=INFOGR.AND.ENTET
C     HTILD : TABLEAU DE TRAVAIL OU ON REFAIT HPROP
C             (TABLEAUX DE MEME STRUCTURE)      
      IF(SOLSYS.EQ.1) THEN
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE CVDFTR SOLSYS=1'
      CALL CVDFTR(T%ADR(ITRAC)%P,TTILD%ADR(ITRAC)%P,TN%ADR(ITRAC)%P,
     *            TSCEXP%ADR(ITRAC)%P,
     *            DIFT,ICONVF(3),CONVV(3),H,HN,HTILD,TETAHC,
     *            UCONV,VCONV,DM1,ZCONV,SOLSYS,
     *            VISCT%ADR(ITRAC)%P,VISC_S,TEXP%ADR(ITRAC)%P,SMH,YASMH,
     *            TIMP%ADR(ITRAC)%P,YASMI(ITRAC),AM1,AM2,ZF,
     *            TBOR%ADR(ITRAC)%P,ATBOR%ADR(ITRAC)%P,
     *            BTBOR%ADR(ITRAC)%P,LIMTRA,MASKTR,MESH,W1,TB,
     *            T1,T2,T3,T4,T5,T6,T7,T10,TE1,TE2,TE3,
     *            KDIR,KDDL,KENT,
     *            DT,ENTET,TETAT,AGGLOT,INFOGT,BILMAS,OPTSUP(3),
     *            ISOUSI,LT,NIT,OPDTRA,OPTBAN,
     *            MSK,MASKEL,MASKPT,MBOR,S,MASSOU(ITRAC),
     *            OPTSOU,SLVTRA,FLBOR,V2DPAR,UNSV2D,2,FLBORTRA)
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE CVDFTR'
C
      ELSE
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE CVDFTR SOLSYS=',SOLSYS
      CALL CVDFTR(T%ADR(ITRAC)%P,TTILD%ADR(ITRAC)%P,TN%ADR(ITRAC)%P,
     *            TSCEXP%ADR(ITRAC)%P,
     *            DIFT,ICONVF(3),CONVV(3),H,HN,HTILD,TETAHC,
     *            UDEL,VDEL,DM1,ZCONV,SOLSYS,
     *            VISCT%ADR(ITRAC)%P,VISC_S,TEXP%ADR(ITRAC)%P,SMH,YASMH,
     *            TIMP%ADR(ITRAC)%P,YASMI(ITRAC),AM1,AM2,ZF,
     *            TBOR%ADR(ITRAC)%P,ATBOR%ADR(ITRAC)%P,
     *            BTBOR%ADR(ITRAC)%P,LIMTRA,MASKTR,MESH,W1,TB,
     *            T1,T2,T3,T4,T5,T6,T7,T10,TE1,TE2,TE3,
     *            KDIR,KDDL,KENT,
     *            DT,ENTET,TETAT,AGGLOT,INFOGT,BILMAS,OPTSUP(3),
     *            ISOUSI,LT,NIT,OPDTRA,OPTBAN,
     *            MSK,MASKEL,MASKPT,MBOR,S,MASSOU(ITRAC),
     *            OPTSOU,SLVTRA,FLBOR,V2DPAR,UNSV2D,2,FLBORTRA)
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE CVDFTR'
      ENDIF
C
      IF(BILMAS) THEN
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE BILANT'
      CALL BILANT(H,T2,T3,DT,LT,NIT,ENTET,
     *            T%ADR(ITRAC)%P,AGGLOT,MASSOU(ITRAC),MASTR0(ITRAC),
     *            MASTR2(ITRAC),MASTEN(ITRAC),
     *            MASTOU(ITRAC),MSK,MASKEL,MESH,FLBOR,NUMLIQ%I,
     *            NFRLIQ,NPTFR,NAMETRAC(ITRAC),FLBORTRA)
      ENDIF
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE BILANT'
C
      ENDDO
C
      ENDIF
C
C=======================================================================
C                   FIN DE LA DIFFUSION DU TRACEUR
C=======================================================================
C
C
C
C
C=======================================================================
C           DIFFUSION ET TERMES SOURCES DU MODELE K-EPSILON
C=======================================================================
C
      IF(ITURB.EQ.3.AND..NOT.ADJO) THEN
C
        IF (ENTET) CALL ENTETE(4,AT,LT)
C
C ATTENTION A LA STRUCTURE DES MATRICES (SYMETRIQUES OU NON)
C QUAND ON PASSERA AU SYSTEME COUPLE K-E
C
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE KEPSIL'
        CALL KEPSIL(AK,EP,AKTILD,EPTILD,AKN,EPN,VISC,CF,U,V,H,
     *              UCONV,VCONV,KBOR,EBOR,LIMKEP%I,IELMK,IELME,
     *              CV1,CV2,TM1,BM1,BM2,CM2,TE1,TE2,NPTFR,DT,
     *              MESH,T1,T2,T3,TB,
     *              CMU,C1,C2,SIGMAK,SIGMAE,ESTAR,SCHMIT,
     *              KMIN,KMAX,EMIN,EMAX,INFOKE.AND.ENTET,
     *              KDIR,MSK,MASKEL,MASKPT,S,SLVK,SLVEP,
     *              ICONVF(4),OPTSUP(4))
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE KEPSIL'
C
      ENDIF
C
C=======================================================================
C  1)                 CONTROLE DU BILAN DE MASSE
C=======================================================================
C
C SECTIONS DE CONTROLE
C
      IF(NCP.NE.0.AND.(ENTET.OR.CUMFLO)) THEN 
        CALL FLUSEC_TELEMAC2D(U,V,H,MESH%IKLE%I,MESH%XEL%R,MESH%YEL%R,
     *                        MESH%NELMAX,MESH%NELEM,
     *                        MESH%X%R,MESH%Y%R,DT,NCP,
     *                        CTRLSC,ENTET,AT,MESH%KNOGL%I,
     *                        MSKSEC,BM1,BM2,T1,HPROP,MESH,S,CV1,
     *                        MESH%IFABOR%I,COMFLU,CUMFLO)
      ENDIF
C
C BILAN DE MASSE
C
      IF(BILMAS) THEN
C
        CALL BILAN(MESH,H,T1,MASK,AT,DT,LT,NIT,ENTET,
     *             MASSES,MSK,MASKEL,EQUA,TE5,OPTBAN,
     *             MESH%NPTFR,FLBOR,
     *             FLUX_BOUNDARIES,NUMLIQ%I,NFRLIQ)
C
C       AJOUTE POUR LES SCHEMAS CINETIQUES (A VERIFIER)
C
        IF(NTRAC.GT.0) THEN
          IF(EQUA(1:15).EQ.'SAINT-VENANT VF') THEN
C
            DO ITRAC=1,NTRAC
            CALL BILANT1(HSTOK,UCONV,VCONV,HPROP,T2,T3,T4,T5,T6,
     *                   DT,LT,NIT,ENTET,MASKTR,
     *                   T%ADR(1)%P,TN%ADR(1)%P,TETAT,
     *                   MASSOU(ITRAC),MSK,MASKEL,MESH,
     *                   FLUTSOR(ITRAC),FLUTENT(ITRAC),EQUA,LTT,ITRAC)
            ENDDO
C
          ENDIF
        ENDIF
C
      ENDIF
C
C=======================================================================
C                        DERIVE DE FLOTTEUR(S)
C=======================================================================
C
      IF(NFLOT.NE.0) THEN
C
        IF(ENTET) CALL ENTETE(12,AT,LT)
C
        IF(SPHERI) THEN
          CALL OS('X=Y/Z   ',UCONV,UCONV,MESH%COSLAT,C)
          CALL OS('X=Y/Z   ',VCONV,VCONV,MESH%COSLAT,C)
        ENDIF
C
        CALL DERIVE(UCONV%R,VCONV%R,DT,
     *              MESH%X%R,MESH%Y%R,MESH%IKLE%I,MESH%IFABOR%I,
     *              LT,IELM,3,NPOIN,
     *              NELEM,NELMAX,MESH%SURDET%R,XFLOT%R,YFLOT%R,
     *              SHPFLO%R,DEBFLO%I,FINFLO%I,ELTFLO%I,
     *              NFLOT,NITFLO,FLOPRD,T8%R)
C
        IF(SPHERI) THEN
          CALL OS('X=XY    ',UCONV,MESH%COSLAT,S,C)
          CALL OS('X=XY    ',VCONV,MESH%COSLAT,S,C)
        ENDIF
C
      ENDIF
C
C=======================================================================
C       OIL SPILL MODEL (UNDER CONSTRUCTION IN MYGRHYCAR PROJECT)
C=======================================================================
C
      IF(SPILL_MODEL) THEN
C
        CALL OIL_SPILL
C
      ENDIF
C
C=======================================================================
C                        DERIVE(S) LAGRANGIENNE(S)
C=======================================================================
C
      IF(NLAG.NE.0) THEN
C
        IF (ENTET) CALL ENTETE(13,AT,LT)
C
        CALL DERLAG(UCONV%R,VCONV%R,DT,MESH%X%R,MESH%Y%R,
     *              MESH%IKLE%I,MESH%IFABOR%I,LT,IELM,3,NPOIN,
     *              NELEM,NELMAX,MESH%SURDET%R,
     *              XLAG%R,YLAG%R,T1%R,T2%R,IT1%I,SHPLAG%R,
     *              DEBLAG%I,FINLAG%I,ELTLAG%I,NLAG,
     *              T7%R,T8%R,MESH%NBOR%I,MESH%NELBOR%I,
     *              MESH%NULONE%I,NPTFR,MSK,MASKEL%R,MASKPT%R,T8%R)
C
      ENDIF
C
C=======================================================================
C                     CONTROLES DE VRAISEMBLANCE
C                   RECHERCHE D'UN ETAT PERMANENT
C=======================================================================
C
      ARRET1=.FALSE.
      IF(VERLIM) THEN
        CALL ISITOK(H%R,H%DIM1,U%R,U%DIM1,V%R,V%DIM1,NTRAC,
     *              T,T%ADR(1)%P%DIM1,
     *              MESH%X%R,MESH%Y%R,BORNES,ARRET1)
C       CORRECTION SUGGESTED BY NOEMIE DURAND (CHC-NRC) 04/01/2006
        IF(NCSIZE.GT.1) THEN
          STOP2=0
          IF(ARRET1) STOP2=1
          STOP2=P_IMAX(STOP2)
          IF(STOP2.EQ.1) ARRET1=.TRUE.
        ENDIF
      ENDIF
      ARRET2=.FALSE.
      IF(STOPER) THEN
        CALL STEADY(H%R,HN%R,H%DIM1,U%R,UN%R,U%DIM1,V%R,VN%R,
     *              V%DIM1,NTRAC,T,TN,T%ADR(1)%P%DIM1,
     *              CRIPER,ARRET2)
C       CORRECTION BY NOEMIE DURAND (CHC-NRC) 04/01/2006
        IF(NCSIZE.GT.1) THEN
          STOP2=0
          IF(ARRET2) STOP2=1
          STOP2=P_IMIN(STOP2)
          ARRET2=.NOT.(STOP2.EQ.0)
        ENDIF
      ENDIF
      IF(ARRET1.OR.ARRET2) THEN
        LEOPRD=1
        LISPRD=1
      ENDIF
C
      ARRET3=.FALSE.
      CALL TRAPSIG()
      IF(BREAKER) ARRET3=.TRUE.
C
      IF(ARRET1.OR.ARRET2.OR.ARRET3) THEN
        LEOPRD=1
        LISPRD=1
      ENDIF
C
C FH-BMD
C=============================================
C     FOR NEW COUPLING
888   CONTINUE
      IF (SISYPHE_CFD) CONSTFLOW_SIS = .TRUE.
999   CONTINUE
C=============================================
C FH-BMD

C
C=======================================================================
C                      IMPRESSION DES RESULTATS
C=======================================================================
C
      IF(ADJO) THEN
C
        IF(T2D_FILES(T2DRBI)%NAME.NE.' '.AND.
     *     INCLU2(ESTIME,'DEBUG')) THEN
          CALL BIEF_DESIMP('SERAFIN ',VARSORA,
     *                     HIST,0,NPOIN,T2D_FILES(T2DRBI)%LU,
     *                     'STD',-AT,LT,LISPRD,1,
     *                     SORLEOA,SORIMPA,MAXVAR,TEXTE,PTINIG,PTINIL)
        ENDIF
C
      ELSE     
C
        IF(CODE(1:7).EQ.'ESTEL3D') THEN
C
C         SAVE THE DEPTH FOR ESTEL3D
          CALL DEPTH_FILL(H%R)
C
C (NOTE THAT OUTPUTS ARE DONE WITHIN ESTEL3D IN COUPLED MODE)
C
        ELSE
          CALL PRERES_TELEMAC2D
          CALL BIEF_DESIMP(T2D_FILES(T2DRES)%FMT,VARSOR,
     *            HIST,0,NPOIN,T2D_FILES(T2DRES)%LU,'STD',AT,LT,
     *            LISPRD,LEOPRD,
     *            SORLEO,SORIMP,MAXVAR,TEXTE,PTINIG,PTINIL)
        ENDIF
C
C
        IF(INCLUS(COUPLING,'DELWAQ')) THEN
C
C         T3 : MODIFIED DEPTH TO TAKE INTO ACCOUNT MASS-LUMPING
C              IN THE CONTINUITY EQUATION
          IF(ABS(1.D0-AGGLOC).GT.1.D-8) THEN
            CALL VECTOR(T3,'=','MASVEC          ',IELMH,
     *                  1.D0-AGGLOC,H ,S,S,S,S,S,MESH,MSK,MASKEL)
            IF(NCSIZE.GT.1) CALL PARCOM(T3,2,MESH)
            CALL OS('X=XY    ',X=T3 ,Y=UNSV2D)
            CALL OS('X=X+CY  ',X=T3 ,Y=H ,C=AGGLOC)
          ELSE
            CALL OS('X=Y     ',X=T3 ,Y=H )
          ENDIF                     
C         FOR COMPUTATION OF THE FLUXES (CALL VECTOR BELOW)
          FORMUL='HUGRADP         '
          IF(SOLSYS.EQ.2) FORMUL(8:8)='2'                  
C           
          CALL VECTOR(T4,'=',FORMUL,11,-1.D0,
     *                HPROP,DM1,ZCONV,UDEL,VDEL,VDEL,MESH,MSK,MASKEL)
          IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE TEL4DEL'
          CALL TEL4DEL(MESH%NPOIN,
     *    MESH%NPOIN,MESH%NELEM,MESH%NSEG,MESH%IKLE%I,MESH%ELTSEG%I,
     *    MESH%GLOSEG%I,MESH%ORISEG%I,MESH%GLOSEG%DIM1,
     *    MESH%X%R,MESH%Y%R,MESH%NPTFR,LIHBOR%I,MESH%NBOR%I,1,
     *    AT,DT,LT,NIT,T3%R,HPROP%R,T3%R,UDEL%R,VDEL%R,
     *    T%ADR(MAX(IND_S,1))%P%R,T%ADR(MAX(IND_T,1))%P%R,
     *    VISC%R,TITCAS,T2D_FILES(T2DGEO)%NAME,
     *    T2D_FILES(T2DCLI)%NAME,WAQPRD,
     *    T2DDL1,T2D_FILES(T2DDL1)%NAME,T2DDL2,T2D_FILES(T2DDL2)%NAME,
     *    T2DDL3,T2D_FILES(T2DDL3)%NAME,T2DDL5,T2D_FILES(T2DDL5)%NAME,
     *    T2DDL6,T2D_FILES(T2DDL6)%NAME,T2DDL7,T2D_FILES(T2DDL7)%NAME,
     *    T2DL11,T2D_FILES(T2DL11)%NAME,T2DDL4,T2D_FILES(T2DDL4)%NAME,
     *    T2DDL8,T2D_FILES(T2DDL8)%NAME,T2DDL9,T2D_FILES(T2DDL9)%NAME,
     *    T2DL10,T2D_FILES(T2DL10)%NAME,ENTET,NELEM,SALI_DEL,
     *    TEMP_DEL,VELO_DEL,DIFF_DEL,MARDAT,MARTIM,FLODEL%R,
     *    INIFLOW,MESH%W%R,.FALSE.,FLULIM%R,V2DPAR%R,MESH%KNOLG%I,
     *    MESH,MESH)
          IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE TEL4DEL'
C
        ENDIF
C
      ENDIF  !(ADJO)
C
C     OPTIONAL USER OUTPUT (COURTESY JACEK JANKOWSKI, BAW)
      IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE UTIMP_TELEMAC2D'
      CALL UTIMP_TELEMAC2D(LT,AT,PTINIG,LEOPRD,PTINIL,LISPRD)
      IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE UTIMP_TELEMAC2D'
C
C=======================================================================
C              COMPARAISON AVEC UN FICHIER DE REFERENCE
C=======================================================================
C
C     LE SOUS-PROGRAMME VALIDA DE LA BIBLIOTHEQUE EST STANDARD.
C     IL PEUT ETRE MODIFIE POUR CHAQUE CAS PARTICULIER.
C     L'APPEL DOIT ETRE LAISSE DANS LA BOUCLE EN TEMPS.
C
      IF(VALID) THEN
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE BIEF_VALIDA'
        CALL BIEF_VALIDA(TB,TEXTPR,
     *                   T2D_FILES(T2DREF)%LU,T2D_FILES(T2DREF)%FMT,
     *                   VARSOR,TEXTE,
     *                   T2D_FILES(T2DRES)%LU,T2D_FILES(T2DRES)%FMT,
     *                   MAXVAR,NPOIN,LT,NIT,ALIRE)
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE BIEF_VALIDA'
      ENDIF
C
C=======================================================================
C
C  CAS OU UN ARRET EST DEMANDE :
C
      IF(ARRET1) THEN
        IF(LNG.EQ.1) THEN
          WRITE(LU,*)
          WRITE(LU,*) 'VALEURS LIMITES DEPASSEES, ARRET DE TELEMAC-2D'
          WRITE(LU,*)
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*)
          WRITE(LU,*) 'LIMIT VALUES TRESPASSED, TELEMAC-2D IS STOPPED'
          WRITE(LU,*)
        ENDIF
        RETURN
      ENDIF
      IF(ARRET2) THEN
        IF(LNG.EQ.1) THEN
          WRITE(LU,*)
          WRITE(LU,*) 'ETAT PERMANENT ATTEINT, ARRET DE TELEMAC-2D'
          WRITE(LU,*)
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*)
          WRITE(LU,*) 'STEADY STATE REACHED, TELEMAC-2D IS STOPPED'
          WRITE(LU,*)
        ENDIF
        RETURN
      ENDIF
      IF(ARRET3) THEN
        IF(LNG.EQ.1) THEN
          WRITE(LU,*)
          CALL ENTETE(1,AT,LT)
          WRITE(LU,*) 'TELEMAC-2D ARRETE PAR L''UTILISATEUR'
          WRITE(LU,*) 'AVEC SIGNAL ',SIGUSR1
          WRITE(LU,*)
        ENDIF
        IF(LNG.EQ.2) THEN
          CALL ENTETE(1,AT,LT)
          WRITE(LU,*)
          WRITE(LU,*) 'TELEMAC-2D CHECKPOINTED BY THE USER'
          WRITE(LU,*) 'USING SIGNAL ',SIGUSR1
          WRITE(LU,*)
        ENDIF
        RETURN
      ENDIF
C
C 700: BOUCLE DES PAS DE TEMPS
C
      IF(LT.LT.NIT) GO TO 700
C
C=======================================================================
C
C :                 /* FIN DE LA BOUCLE EN TEMPS */
C
C=======================================================================
C
      IF(LNG.EQ.1.AND.LISTIN) WRITE(LU,18)
      IF(LNG.EQ.2.AND.LISTIN) WRITE(LU,19)
18    FORMAT(/,1X,'FIN DE LA BOUCLE EN TEMPS',////)
19    FORMAT(/,1X,'END OF TIME LOOP',////)
C
C-----------------------------------------------------------------------
C
      IF (NFLOT.NE.0) CALL SORFLO
     *   (XFLOT%R,YFLOT%R,IKLFLO%I,DEBFLO%I,FINFLO%I,
     *    NFLOT,NITFLO,FLOPRD,T2D_FILES(T2DRBI)%LU,TITCAS,
     *    'STD',T2D_FILES(T2DRBI)%NAME,
     *    NIT,MAXVAR,MARDAT,MARTIM,MESH,I_ORIG,J_ORIG)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
