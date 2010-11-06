C                       ******************
                        SUBROUTINE SISYPHE
C                       ******************
C
     *(PART,LOOPCOUNT,GRAFCOUNT,LISTCOUNT,TELNIT,
     * U_TEL,V_TEL,H_TEL,HN_TEL,ZF_SIS,UETCAR,CF_TEL,
     * CONSTFLOW,NSIS_CFD,SISYPHE_CFD,CODE,PERICOU,
     * U3D,V3D,T_TEL,VISC_TEL,DT_TEL,CHARR_TEL,SUSP_TEL,
     * FLBOR_TEL,SOLSYS,DM1,UCONV_TEL,VCONV_TEL,ZCONV)
C
C**********************************************************************
C SISYPHE VERSION 6.0                             C. LENORMANT   
C                                                 J.-M. HERVOUET
C                                                 S. HADJI
C                                                 C. MACHET
C                                                 C. VILLARET
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
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____._______________________________________________
C |      NOM       |MODE|                   ROLE                        
C |________________|____|_______________________________________________
C |   PART         | -->| SI -1, PAS DE COUPLAGE : ON PASSE TOUTE LA
C |                |    | SUBROUTINE. SINON, INDIQUE LA PARTIE DE LA
C |                |    | SUBROUTINE DANS LAQUELLE ON PASSE
C |                |    |
C  LES ARGUMENTS SUIVANTS NE SERVENT QUE EN COUPLAGE AVEC TELEMAC:
C |   LOOPCOUNT    | -->| NUMERO DE L'ITERATION        
C |   GRAFCOUNT    | -->| PERIODE DE SORTIE GRAPHIQUE    
C |   LISTCOUNT    | -->| PERIODE DE SORTIE LISTING                               
C |   TELNIT       | -->| NOMBRE D'ITERATION                                              
C |   U_TEL        | -->|
C |   V_TEL        | -->| VARIABLES HYDRO ENVOYEES PAR TELEMAC 2D 
C |   H_TEL        | -->|
C |   ZF_SIS       |<-- | FOND ENVOYE A TELEMAC 2D   
C |   UETCAR       | -->| VARIABLES HYDRO ENVOYEES PAR TELEMAC 2D 
C |   CONSTFLOW    |    |
C |   NSIS_CFD     |    |
C |   SISYPHE_CFD  |    |
C |   CODE         |    | NAME OF CALLING PROGRAMME
C |   PERICOU      |    |
C |   U3D,V3D      | -->| 3D VELOCITY SENT BY TELEMAC 3D 
C |   T_TEL        | -->| CURRENT TIME IN CALLING PROGRAMME
C----------------------------------------------------------------------
C PROGRAMME APPELANT : HOMERE
C PROGRAMMES APPELES : PREDAT , ENTETE , CHECK  , LECT   ,
C                      CONDIM , PREDES , DESIMP , CORSTR , TRANSP ,
C                      TOPOGR ,ARRET  , MASKBD , MASKAB  , CFL    ,
C                      MAXI   , MINI   , CONLIT , RESOLU , VALIDA ,
C                      SP DE BIEF. + FINREF , DESREF, CALCUL_SUSP
C
C
C
!======================================================================!
!======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS               !
!======================================================================!
!======================================================================!
!
      USE INTERFACE_SISYPHE, EX_SISYPHE => SISYPHE
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_SISYPHE
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER,           INTENT(IN)    :: PART,LOOPCOUNT,GRAFCOUNT
      INTEGER,           INTENT(IN)    :: LISTCOUNT,TELNIT,PERICOU
      CHARACTER(LEN=24), INTENT(IN)    :: CODE
      TYPE(BIEF_OBJ),    INTENT(IN)    :: U_TEL,V_TEL,H_TEL,HN_TEL
      TYPE(BIEF_OBJ),    INTENT(INOUT) :: ZF_SIS,UETCAR     
      INTEGER,           INTENT(INOUT) :: NSIS_CFD
      LOGICAL,           INTENT(INOUT) :: CONSTFLOW,SISYPHE_CFD 
      TYPE(BIEF_OBJ),    INTENT(IN)    :: U3D,V3D,VISC_TEL
      TYPE(BIEF_OBJ),    INTENT(INOUT) :: CF_TEL
      DOUBLE PRECISION,  INTENT(IN)    :: T_TEL
      LOGICAL,           INTENT(INOUT) :: CHARR_TEL,SUSP_TEL
      DOUBLE PRECISION,  INTENT(IN)    :: DT_TEL
      INTEGER,           INTENT(IN)    :: SOLSYS
      TYPE(BIEF_OBJ),    INTENT(IN)    :: FLBOR_TEL,DM1,ZCONV
      TYPE(BIEF_OBJ),    INTENT(IN)    :: UCONV_TEL,VCONV_TEL
C
      INTEGER                        P_IMAX
      DOUBLE PRECISION P_DMAX,P_DMIN
      EXTERNAL         P_DMAX,P_DMIN,P_IMAX
C
      INTEGER, PARAMETER :: NHIST = 0
      INTEGER, PARAMETER :: NSOR = 100
      INTEGER            :: VALNIT,NLISS
      INTEGER            :: I,J,K,MN,MT,ISOUS,NIDT,NIT,IMA,IMI
      INTEGER            :: IMIN,IMAX,NCALCU,NUMEN,NUMENX,NUMDEB
      INTEGER            :: ALIRE(MAXVAR),ALIRV(MAXVAR),ALIRH(MAXVAR)
      INTEGER            :: ALIR0(MAXVAR)
      INTEGER            :: TROUVE(MAXVAR+10)
      DOUBLE PRECISION   :: AT0,DTS,BID,XMA,XMI
      DOUBLE PRECISION   :: XMIN,XMAX
      DOUBLE PRECISION   :: AT,VCUMU,MASS_GF
      DOUBLE PRECISION   :: HIST(1)
      LOGICAL            :: PASS,PASS_SUSP
      LOGICAL            :: RESU,ENTETS,CHGMSK,YAZR
!
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: SAVEZF,SAVEQU,SAVEQV
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: SAVEZ
!
      ! SAUVEGARDE DES VARIABLES LOCALES
      ! --------------------------------
      SAVE VCUMU               ! POUR LE BILAN
      SAVE MASS_GF             ! FOR GRAIN-FEEDING
      SAVE PASS, PASS_SUSP     ! MARQUE DU PREMIER PAS DE TEMPS
      SAVE NIDT, NCALCU, NUMEN, NIT, VALNIT ! 
      SAVE AT0                 ! TEMPS
C     NUMEN0 : PREMIER ENREGISTREMENT A LIRE
      INTEGER :: NUMEN0
!
      ! VARIABLES A LIRE EN CAS DE SUITE
      ! --------------------------------
      ! 0 : ON NE LIT PAS
      ! 1 : ON LIT  (VOIR SS-PG NOMVAR)
!
C   HYDRO + EVOLUTION
      DATA ALIRE /1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     *            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     *            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     *            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
! HOULE SEULE
      DATA ALIRH /0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
     *            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     *            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     *            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
! NE RIEN LIRE
      DATA ALIR0 /0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     *            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     *            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     *            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
!
! FOR VALIDATION, EVERY VARIABLE IN THE FILE IS COMPARED
!
      DATA ALIRV /1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     *            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     *            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     *            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/
!
!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!
!
!
C------------------------------------------------------------------
C     PARTIE 1 : INITIALIZATION
C------------------------------------------------------------------
C    
      IF(PART==0.OR.PART==-1) THEN
        IF(DEBUG.GT.0) WRITE(LU,*) 'INITIALIZATION' 
C
        WRITE(LU,*) 'PART 0 : INITIALISING SISYPHE'
C        
C       INITIALISATION DES CONSTANTES
C       
        CALL INIT_CONSTANT(KARIM_HOLLY_YANG,KARMAN,PI)
C        IF(SUSP) CALL USER_KRONE_PART(VITCE,VITCD,
C     *                                PARTHENIADES,XMVS,CMAX)
C          
C      LECTURE DES DONNEES DE HOULE SUR FICHIER HYDRODYNAMIQUE
C
        IF(HOULE.AND.SIS_FILES(SISCOU)%NAME(1:1).EQ.' ') THEN
          ALIRE(12)=1
          ALIRE(13)=1
          ALIRE(14)=1
        ENDIF
C
C       LECTURE DES DONNEES SEDIMENTO SUR FICHIER SUITE
C
        IF(DEBU) THEN
          ALIRE(15)=1
          ALIRE(16)=1
          ALIRE(17)=1
          ALIRE(18)=1                       
C         TO READ THE AVAI FROM THE PREVIOUS COMPUTATION FILE
          DO I=1,NSICLA*NOMBLAY
            ALIRE(20+I)=1
          ENDDO           
C         TO READ THE CS (CONCENTRATION) FROM THE PREVIOUS COMPUTATION FILE
          IF(SUSP) THEN
            DO I=1,NSICLA
             ALIRE(20+(NOMBLAY+1)*NSICLA+I)=1
            ENDDO 
          ENDIF
C         TO READ THE LAYER THICKNESSES
          DO I=1,NOMBLAY
            ALIRE(26+(NOMBLAY+4)*NSICLA+I)=1
          ENDDO
        ENDIF
C
C --------  MISE A ZERO DES TABLEAUX
C
        CALL INIT_ZERO   
C
C --------  FIN MISE A ZERO DES TABLEAUX
C
C       DISCRETISATIONS : LINEAIRE POUR L'INSTANT.
C
C       IELMT CODE EN DUR DANS LECDON
C
        IELMX  = MAX(IELMT,IELMH_SIS,IELMU_SIS)
        NELMAX = NELEM
C
C=======================================================================
C
C : 1          LECTURE , PREPARATION ET CONTROLE DES DONNEES
C
C=======================================================================
C
        RESU = .TRUE.
C
C       LECTURE DES CONDITIONS LIMITES ET INDICES DES POINTS FRONTIERES.
C
        IF(DEBUG.GT.0) WRITE(LU,*) 'LECLIS'
        CALL LECLIS(LIEBOR%I,EBOR,
     &              MESH%NPTFR,MESH%NBOR%I,3,
     &              SIS_FILES(SISCLI)%LU,KENT,
     &              MESH%ISEG%I,MESH%XSEG%R,MESH%YSEG%R,MESH%NACHB%I,
     &              NUMLIQ%I,NSICLA,AFBOR%R,BFBOR%R,BOUNDARY_COLOUR%I,
     &              MESH)
        IF(DEBUG.GT.0) WRITE(LU,*) 'END_LECLIS'
C     
C-----------------------------------------------------------------------
C
C       COMPLEMENT DE LA STRUCTURE DE DONNEES POUR BIEF
C
        IF(DEBUG.GT.0) WRITE(LU,*) 'INBIEF'
        CALL INBIEF(IT1%I,KLOG,IT2,IT3,IT4,LVMAC,IELMX,
     *                 0.D0,SPHERI,MESH,T1,T2,OPTASS,PRODUC,EQUA)
C
        IF(DEBUG.GT.0) WRITE(LU,*) 'END_INBIEF'
C
C-----------------------------------------------------------------------
C
C       REPERAGE DES FRONTIERES
C
        IF(NCSIZE.GT.1) THEN
          NFRLIQ=0
          DO I=1,NPTFR
            NFRLIQ=MAX(NFRLIQ,NUMLIQ%I(I))
          ENDDO
          NFRLIQ=P_IMAX(NFRLIQ)
          WRITE(LU,*) ' '
          IF(LNG.EQ.1) WRITE(LU,*) 'FRONTIERES LIQUIDES :',NFRLIQ
          IF(LNG.EQ.2) WRITE(LU,*) 'LIQUID BOUNDARIES:',NFRLIQ
        ELSE
          IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE FRONT2'
          CALL FRONT2(NFRLIQ,NFRSOL,DEBLIQ,FINLIQ,DEBSOL,FINSOL,
     *                LIEBOR%I,LIEBOR%I,
     *                MESH%X%R,MESH%Y%R,MESH%NBOR%I,MESH%KP1BOR%I,
     *                IT1%I,NPOIN,NPTFR,KLOG,.TRUE.,NUMLIQ%I,MAXFRO)
          IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE FRONT2'
        ENDIF
C
C-----------------------------------------------------------------------
C       RECHERCHE DU FOND ET DU FROTTEMENT DANS LE FICHIER DE GEOMETRIE:
C-----------------------------------------------------------------------
C
        IF(DEBUG.GT.0) WRITE(LU,*) 'FONSTR'
        CALL FONSTR(T1,ZF,T2,CHESTR,SIS_FILES(SISGEO)%LU,
     *              SIS_FILES(SISFON)%LU,SIS_FILES(SISFON)%NAME,
     *              MESH,SFON,.TRUE.)
        IF(DEBUG.GT.0) WRITE(LU,*) 'END_FONSTR'
C
C-----------------------------------------------------------------------
C
C       PREPARATION DU FICHIER DE RESULTATS (FACULTATIF)
C
C       STANDARD SELAFIN
C
        IF(DEBUG.GT.0) WRITE(LU,*) 'ECRGEO'
        ! CREATION DU JEU DE DONNEES POUR UN FORMAT DE FICHIER
        ! FORMAT_RES.
        ! LE JEU DE DONNEES EST CREE DANS LE FICHIER SISRES, ET EST
        ! DEFINIT PAR UN TITRE ET LES VARIABLES A ECRIRE. 
        CALL CREATE_DATASET(SIS_FILES(SISRES)%FMT, ! FORMAT FICHIER RESULTAT
     *                      SIS_FILES(SISRES)%LU,  ! LU FICHIER RESULTAT
     *                      TITCA,      ! TITRE DE L'ETUDE
     *                      MAXVAR,     ! MAX VARIABLES SORTIE
     *                      TEXTE,      ! NOMS VARIABLES SORTIE
     *                      SORLEO)     ! SORTIE OU PAS DES VARIABLES
        ! ECRITURE DU MAILLAGE DANS LE FICHIER SORTIE :
        ! SI ON EST ON PARALLEL, FAUT L'INDIQUER VIA NCSIZE ET NPTIR.
        ! LES AUTRES INFORMATIONS SONT DANS MESH.
        ! EN PLUS : DATE/TEMPS DE DEPART ET LES COORDONNEES DE
        ! L'ORIGINE.
        CALL WRITE_MESH(SIS_FILES(SISRES)%FMT, ! FORMAT FICHIER RESULTAT     
     *                  SIS_FILES(SISRES)%LU,  ! LU FICHIER RESULTAT
     *                  MESH,          ! DESCRIPTEUR MAILLAGE
     *                  1,             ! NOMBRE DE PLAN /NA/
     *                  MARDAT,        ! DATE DEBUT
     *                  MARTIM,        ! HEURE DEBUT
     *                  I_ORIG,J_ORIG) ! COORDONNEES DE L'ORIGINE.
        IF(DEBUG.GT.0) WRITE(LU,*) 'END_ECRGEO'
C
C   --- REMPLISSAGE PAR DEFAUT DU TABLEAU MASKEL
C
        IF(MSK) CALL OS ('X=C     ', X=MASKEL, C=1.D0)
C
C       CONSTRUCTION DU MASQUE
C
        DO K = 1, MESH%NPTFR
          IF(LIEBOR%I(K).NE.2.AND.LIEBOR%I(MESH%KP1BOR%I(K)).NE.2) THEN
            MASK%R(K) = 1.D0
          ELSE
            MASK%R(K) = 0.D0
          ENDIF
          LIQBOR%I(K) = KSORT
        ENDDO
C
C=======================================================================
C
C : 2                  INITIALISATIONS
C
C=======================================================================
C
        PASS      = .TRUE.
C        
        PASS_SUSP = .TRUE.
        VCUMU     = 0.D0
        MASS_GF   = 0.D0
C
C 
C----   DETERMINATION DU NOMBRE D'EVENEMENTS (1ere boucle)       : NCALCU 
C                               NOMBRE DE PAS DE TEMPS (2e boucle): NIDT 
C                               NOMBRE TOTAL DE PAS DE TEMPS     : NIT 
c        	                TEMPS INITIAL                    : AT0
C                    		PAS DE TEMPS DE CALCUL           : DT 
C
C
        IF(PART.EQ.0) THEN
           AT0=T_TEL
           DT = DT_TEL
           NCALCU = 1
           NIDT   = 1
           NIT=TELNIT
        ELSE
          AT0=MAX(TPREC,0.D0)
          DT=DELT
          IF(PERMA) THEN 
             NCALCU=1
             NIDT=NPAS
             NIT=NIDT
          ELSE 
             NCALCU = NMAREE 
C DT calcule apres lecture fichier hydro
C             NIDT =  NINT ( PMAREE / DT + 0.1D0 )
C             NIT=NIDT*NCALCU
C sinon 
              NIDT=NPAS
              NIT=NIDT*NCALCU
          ENDIF
        ENDIF  
C
C En nom permanent : DT est calcule a partir du Fichier Hydro
C                    NUMEN: nombre total d'enregistrement

        IF(SIS_FILES(SISHYD)%NAME(1:1).NE.' ')  THEN
            IF(DEBUG.GT.0) WRITE(LU,*) 'BIEF_SUITE'
            WRITE(LU,*) 'APPEL DE BIEF_SUITE'
!           JUST TO GET NUMEN AND DT (SEE ALIR0)
            CALL BIEF_SUITE(VARSOR,VARCL,NUMEN,SIS_FILES(SISHYD)%LU,
     *                    SIS_FILES(SISHYD)%FMT,HIST,0,NPOIN,AT,
     *                    TEXTPR,VARCLA,
     *                    0,TROUVE,ALIR0,.TRUE.,.TRUE.,MAXVAR,DT=DT)
            IF(DEBUG.GT.0) WRITE(LU,*) 'SORTIE DE BIEF_SUITE'
            WRITE(LU,*) 'LECTURE FICHIER HYDRODYNAMIQUE:'
          IF(.NOT.PERMA) THEN
             NIDT =  NINT ( PMAREE / DT + 0.1D0 )
             IF(ABS(NIDT*DT-PMAREE) > 1.D-3) THEN
               IF (LNG == 1) WRITE(LU,101) NIDT*DT
               IF (LNG == 2) WRITE(LU,102) NIDT*DT
             ENDIF   
             NIT  = NCALCU * NIDT
          ENDIF          
        ENDIF  
C
C Validation au dernier pas de temps
C                
        VALNIT = NIT
C
101     FORMAT(/,'ATTENTION : LA PERIODE DE CALCUL NE CORRESPOND PAS A',/,
     *       'UN MULTIPLE DE LA PERIODE DE SORTIE HYDRODYNAMIQUE.',/,
     *       'LE CALCUL S''EFFECTUERA DONC SUR ',G16.7,'SECONDES')
102     FORMAT(/,'CAUTION : THE PERIOD OF COMPUTATION IS NOT A MULTIPLE',
     *         /,'OF THE HYDRODYNAMIC FILE PRINTOUT PERIOD.',/,
     *       'THE LENGTH OF COMPUTATION WILL THEREFORE BE',G16.7,/,
     *       'SECONDS')     
C

C  MODE SISYPHE SEUL
C  -----------------------------------------------------------------------
C  ---- LECTURE 
C       DU FICHIER PRECEDENT HYDRODYNAMIQUE
C 
C
C NUMEN : NBRE D'ENREGISTREMENTS DU FICHIER HYDRODYNAMIQUE
C DT :  DETERMINER PAS DE TEMPS ENREGISTREMENTS HYDRODYNAMIQUES
C NUMEN0= 1ER ENREGISTREMENT A LIRE SUR FICHIER HYDRO
C TPREC : TEMPS DE DEMARRAGE
C
C 
C V5P9      NUMEN0 = INT( (TPREC - ATDEB)/DT + 1.1D0 )
C              
        IF(PART.EQ.-1) THEN
          IF(.NOT.PERMA) THEN
            IF(TPREC.GE.0.D0) THEN
                NUMEN0 = INT( TPREC /DT + 1.1D0 )
            ELSE   
                NUMEN0 = NUMEN-INT(PMAREE/DT+1.1D0)
            ENDIF    
          ELSE  
            IF(TPREC.GE.0.D0) THEN
               NUMEN0 = INT( TPREC /DT + 1.1D0 )
            ELSE   
               NUMEN0 = NUMEN
            ENDIF    
          ENDIF
C
          IF(SIS_FILES(SISHYD)%NAME(1:1).NE.' ')  THEN
            IF(DEBUG.GT.0) WRITE(LU,*) 'BIEF_SUITE'
            CALL BIEF_SUITE(VARSOR,VARCL,NUMEN0,SIS_FILES(SISHYD)%LU,
     *                    SIS_FILES(SISHYD)%FMT,HIST,0,NPOIN,AT,
     *                    TEXTPR,VARCLA,
     *                    0,TROUVE,ALIRE,.TRUE.,PERMA,MAXVAR)
C
C         TRACING IF WAVE DATA HAVE BEEN FOUND
C
            IF(HOULE) THEN
              IF(TROUVE(12).EQ.1) HW%TYPR='Q'
              IF(TROUVE(13).EQ.1) TW%TYPR='Q'
              IF(TROUVE(14).EQ.1) THETAW%TYPR='Q'
            ENDIF
            IF(DEBUG.GT.0) WRITE(LU,*) 'END_BIEF_SUITE'
            IF(DEBUG.GT.0) WRITE(LU,*) 'RESCUE_SISYPHE'
            CALL RESCUE_SISYPHE(QU%R,QV%R,Q%R,U2D%R,V2D%R,HN%R,Z%R,
     *                        ZF%R,HW%R,TW%R,THETAW%R,NPOIN,
     *                        TROUVE,ALIRE,PASS,ICF,.TRUE.)
            IF(DEBUG.GT.0) WRITE(LU,*) 'END_RESCUE_SISYPHE'
          ENDIF

C
        ENDIF
C 
C---- SUITE D'UN CALCUL SISYPHE
C
        YAZR=.FALSE.     
        IF(SIS_FILES(SISPRE)%NAME(1:1).NE.' ')  THEN      
C
C         LECTURE DES VARIABLES HYDRO ET SEDIMENTOLOGIQUES
C
          IF(DEBUG.GT.0) WRITE(LU,*) 'BIEF_SUITE'
          CALL BIEF_SUITE(VARSOR,VARCL,NUMENX,SIS_FILES(SISPRE)%LU,
     *                    SIS_FILES(SISPRE)%FMT,   
     *                    HIST,0,NPOIN,AT0,TEXTPR,VARCLA,0,
     *                    TROUVE,ALIRE,.TRUE.,.TRUE.,MAXVAR)
          IF(TROUVE(9).EQ.1) YAZR=.TRUE.
          IF(DEBUG.GT.0) WRITE(LU,*) 'END_BIEF_SUITE'
C
          IF(DEBUG.GT.0) WRITE(LU,*) 'RESCUE_SISYPHE'
          CALL RESCUE_SISYPHE(QU%R,QV%R,Q%R,U2D%R,V2D%R,HN%R,Z%R,ZF%R,
     *                        HW%R,TW%R,THETAW%R,NPOIN,TROUVE,ALIRE,
     *                        PASS,ICF,.TRUE.)
          IF(DEBUG.GT.0) WRITE(LU,*) 'SORTIE DE BIEF_SUITE'
C
C         CHANGING UNITS OF CONCENTRATIONS
C
          IF(SUSP.AND.UNIT) THEN
            DO I=1,NSICLA
             IF(TROUVE(20+(NOMBLAY+1)*NSICLA+I).EQ.1) THEN
               CALL OS('X=CX    ',X=CS%ADR(I)%P,C=1.D0/XMVS)
             ENDIF
            ENDDO
          ENDIF
C
        ENDIF
C
C----    LECTURE DU DERNIER ENEGISTREMENT : FICHIER DE HOULE
C      
C       NOTE : SIS_FILES(SISCOU)%NAME SET TO ' ' IF HOULE=NO
C
        IF(SIS_FILES(SISCOU)%NAME(1:1).NE.' ')  THEN
C
          IF(DEBUG.GT.0) WRITE(LU,*) 'SUITE_HOULE'
          WRITE(LU,*) ' LECTURE HOULE :',SIS_FILES(SISCOU)%NAME
          CALL BIEF_SUITE(VARSOR,VARCL,NUMENX,SIS_FILES(SISCOU)%LU,
     *                    SIS_FILES(SISCOU)%FMT,HIST,0,NPOIN,AT,
     *                    TEXTPR,VARCLA,0,
     *                    TROUVE,ALIRH,.TRUE.,.TRUE.,MAXVAR)
          IF(DEBUG.GT.0) WRITE(LU,*) 'END_SUITE_HOULE'
C         TRACING IF WAVE DATA HAVE BEEN FOUND
          IF(TROUVE(12).EQ.1) HW%TYPR='Q'
          IF(TROUVE(13).EQ.1) TW%TYPR='Q'
          IF(TROUVE(14).EQ.1) THETAW%TYPR='Q'
C 
        ENDIF
C 
        IF(CODE(1:7) == 'TELEMAC'.AND.PART==0) THEN
          WRITE(LU,*) 'INITIALISATION EN CAS DE COUPLAGE : PART=',PART
C         INFORMATION ON SUSPENSION SENT BACK
          CHARR_TEL = CHARR 
          SUSP_TEL = SUSP
C
C         OV INSTEAD OF OS IN ORDER TO AVOID PROBLEMS WITH QUASI-BUBBLE ELEMENTS
C         ONE OPERATES ONLY ON THE (1:NPOIN) RANGE OF THE TELEMAC FIELDS
C         IT IS A HIDDEN DISCRETISATION CHANGE
C
          CALL OV( 'X=Y     ', U2D%R, U_TEL%R, U_TEL%R, 0.D0, NPOIN)
          CALL OV( 'X=Y     ', V2D%R, V_TEL%R, V_TEL%R, 0.D0, NPOIN)
          CALL OV( 'X=Y     ', HN%R , H_TEL%R, H_TEL%R, 0.D0, NPOIN)
          CALL OS('X=Y     ', X=ZF,Y=ZF_SIS)
C
C         CLIPPING OF NEGATIVE DEPTHS
C
          IF(OPTBAN.GT.0) THEN
            DO I = 1,NPOIN
             IF(HN%R(I).LT.HMIN) THEN
               U2D%R(I)=0.D0
               V2D%R(I)=0.D0
               HN%R(I) = HMIN
             ENDIF            
            ENDDO
          ENDIF
        ENDIF
C
C  ---- END COUPLING  -------------
C
        IF(DEBUG.GT.0) WRITE(LU,*) 'CONDIM_SISYPHE'
        IF(.NOT.DEBU) THEN
        CALL CONDIM_SISYPHE
     *        (U2D%R,V2D%R,QU%R,QV%R,HN%R,ZF%R,Z%R,ESOMT%R,THETAW%R,
     *         Q%R,HW%R,TW%R,MESH%X%R,MESH%Y%R,NPOIN,AT0,PMAREE)
        ENDIF
        IF(DEBUG.GT.0) WRITE(LU,*) 'END_CONDIM_SISYPHE'                
C
C       AT THIS LEVEL U2D,V2D,H AND ZF MUST HAVE BEEN DEFINED
C       EITHER BY BIEF_SUITE, CONDIM_SISYPHE OR CALLING PROGRAMME    
C
C       NOW COMPUTING FUNCTIONS OF U2D,V2D,HN AND ZF
C
C       FREE SURFACE
        CALL OS('X=Y+Z   ', X=Z, Y=ZF, Z=HN) 
C
        IF(CODE(1:7).NE.'TELEMAC') THEN       
C         PRODUCT H*
          CALL OS('X=YZ    ', X=QU, Y=U2D, Z=HN)
C         PRODUCT H*V
          CALL OS('X=YZ    ', X=QV, Y=V2D, Z=HN)
C         DISCHARGE 
          CALL OS('X=N(Y,Z)', X=Q, Y=QU, Z=QV)
        ENDIF
C
C       CHECKING WAVE DATA
C 
        IF(HOULE) THEN
          IF(HW%TYPR    .NE.'Q'.OR.
     *       TW%TYPR    .NE.'Q'.OR.
     *       THETAW%TYPR.NE.'Q') THEN
            WRITE(LU,*) ' '
            WRITE(LU,*) ' '
            IF(LNG.EQ.1) THEN
              WRITE(LU,*) 'DONNEES DE HOULE MANQUANTES'
              IF(HW%TYPR.NE.'Q') WRITE(LU,*) 'HAUTEUR HM0'
              IF(TW%TYPR.NE.'Q') WRITE(LU,*) 'PERIODE PIC TPR5'
              IF(THETAW%TYPR.NE.'Q') WRITE(LU,*) 'DIRECTION MOY'
            ENDIF 
            IF(LNG.EQ.2) THEN
              WRITE(LU,*) 'MISSING WAVE DATA'
              IF(HW%TYPR.NE.'Q') WRITE(LU,*) 'WAVE HEIGHT HM0'
              IF(TW%TYPR.NE.'Q') WRITE(LU,*) 'PEAK PERIOD TPR5'
              IF(THETAW%TYPR.NE.'Q') WRITE(LU,*) 'MEAN DIRECTION'
            ENDIF
            CALL PLANTE(1)
            STOP
          ENDIF
        ENDIF
C
C FIN INITIALISATION HYDRODYNAMIQUE
C
C
C        CALCUL DES SURFACES (SANS MASQUAGE) 
C
         CALL VECTOR(VOLU2D,'=','MASBAS          ',
     *               IELMH_SIS,1.D0, 
     *               T1,T1,T1,T1,T1,T1,MESH,.FALSE.,MASKEL)
C        V2DPAR : LIKE VOLU2D BUT IN PARALLEL VALUES COMPLETED AT 
C                 INTERFACES BETWEEN SUBDOMAINS
         CALL OS('X=Y     ',X=V2DPAR,Y=VOLU2D)
         IF(NCSIZE.GT.1) CALL PARCOM(V2DPAR,2,MESH)
C        INVERSE OF VOLUMES (DONE WITHOUT MASKING)
         CALL OS('X=1/Y   ',X=UNSV2D,Y=V2DPAR,
     *           IOPT=2,INFINI=0.D0,ZERO=1.D-12)
C
C DEBUT MODIFS SEDIMENTS MIXTES
C
C        SETTING THE NON ERODABLE BED (IT CAN BE SET BEFORE
C                                      IF COMPUTATION CONTINUED, I.E. DEBU)
C
         IF(.NOT.DEBU.OR..NOT.YAZR) THEN
           CALL NOEROD(HN%R,ZF%R,ZR%R,Z%R,MESH%X%R,
     *                 MESH%Y%R,NPOIN,CHOIX,NLISS)
         ENDIF    
C
C        INITIALISATION POUR LE SEDIMENT
C
        IF(DEBUG.GT.0) WRITE(LU,*) 'INIT_SEDIMENT'
          CALL INIT_SEDIMENT(NSICLA,ELAY,ZF,ZR,NPOIN,
     *                    AVAIL,FRACSED_GF,AVA0,LGRAFED,CALWC,
     *                    XMVS,XMVE,GRAV,VCE,XWC,FDM,CALAC,AC,
     *                    SEDCO, ES, NCOUCH_TASS,CONC_VASE,
     *                    MS_SABLE, MS_VASE,ACLADM, UNLADM)
          IF(DEBUG.GT.0) WRITE(LU,*) 'END INIT_SEDIMENT' 
C
C
C FIN MODIFS CV 
C
! 
! VITESSE MOYENNE
!======================================================================

        CALL OS('X=N(Y,Z)',X=UNORM,Y=U2D,Z=V2D)         

! =====================================================================
!  VITESSE ORBITALE SOUS LA HOULE
! =====================================================================

        IF(HOULE) THEN
          CALL CALCUW(UW%R,HN%R,HW%R,TW%R,GRAV,NPOIN)
        ENDIF
        
C ======================================================================
        IF(DEBUG.GT.0) WRITE(LU,*) 'TOB_SISYPHE'
        CALL TOB_SISYPHE(TOB,TOBW,MU,KS,KSP,KSR,CF,FW, 
     &                   CHESTR,UETCAR,CF_TEL,CODE ,
     &                   KFROT,ICR,KSPRATIO,HOULE,
     &                   GRAV,XMVE,XMVS,VCE,KARMAN,ZERO,
     &                   HMIN,HN,ACLADM,UNORM,UW,TW,NPOIN)
        IF(DEBUG.GT.0) WRITE(LU,*) 'END TOB_SISYPHE'                                
C
C       INITIALISATION POUR LE TRANSPORT
C
        IF(DEBUG.GT.0) WRITE(LU,*) 'INIT_TRANSPORT'
        CALL INIT_TRANSPORT(TROUVE,DEBU,HIDING,NSICLA,NPOIN,
     *     T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T14,
     *     CHARR,QS_C,QSXC,QSYC,CALFA,SALFA,COEFPN,SLOPEFF,
     *     SUSP,QS_S,QS,QSCL,QSCL_C,QSCL_S,QSCLXS,QSCLYS, 
     *     UNORM,U2D,V2D,HN,CF,MU,TOB,TOBW,UW,TW,THETAW, FW,HOULE,
     *     AVAIL,ACLADM,UNLADM,KSP,KSR,
     *     ICF,HIDFAC,XMVS,XMVE,GRAV,VCE,XKV,HMIN, KARMAN,
     *     ZERO,PI,AC,IMP_INFLOW_C,ZREF,ICQ,CSTAEQ,
     *     CMAX,CS,CS0,UCONV,VCONV,CORR_CONV,SECCURRENT,BIJK,
     *     IELMT, MESH, FDM,XWC,FD90,SEDCO,VITCE,PARTHENIADES,VITCD)
        IF(DEBUG.GT.0) WRITE(LU,*) 'END INIT_TRANSPORT'
C
C ---------- DEBUT IMPRESSIION INITIALISATION =================
C 
        CALL ENTETE_SISYPHE(1,AT0,0)
C       PREPARATION DES RESULTATS         
C
C CONCENTRATION OUTPUT IN G/L
C
        IF(UNIT) CALL OS('X=CX    ',X=CS,C=XMVS)      
        CALL PREDES(0,AT0)
C        
C       IMPRESSION DES RESULTATS
C
        IF(DEBUG.GT.0) WRITE(LU,*) 'BIEF_DESIMP'
        CALL BIEF_DESIMP(SIS_FILES(SISRES)%FMT,VARSOR,
     *                   HIST,0,NPOIN,SIS_FILES(SISRES)%LU,'STD',
     *                   AT0,0,LISPR,LEOPR,SORLEO,SORIMP,MAXVAR,
     *                   TEXTE,PTINIG,PTINIL)
        IF(DEBUG.GT.0) WRITE(LU,*) 'END BIEF_DESIMP'
C
        IF(UNIT) CALL OS('X=CX    ',X=CS,C=1.D0/XMVS)      
C
C===============FIN IMPRESSION CONDITIONS INITIALES =================
C
C      COUPLING
C
       IF(DREDGESIM) THEN
         CALL DREDGESIM_INTERFACE(1)
         IF(LNG.EQ.1) WRITE(LU,*) 'SISYPHE COUPLE AVEC DREDGESIM'
         IF(LNG.EQ.2) WRITE(LU,*) 'SISYPHE COUPLED WITH DREDGESIM'
       ENDIF
C
       IF(CODE(1:7).NE.'SISYPHE') THEN
         IF(LNG.EQ.1) WRITE(LU,*) 'SISYPHE COUPLE AVEC : ',CODE
         IF(LNG.EQ.2) WRITE(LU,*) 'SISYPHE COUPLED WITH: ',CODE
       ENDIF
C
C      COUPLING WITH TELEMAC-2D OR 3D
C
       IF(CODE(1:7).EQ.'TELEMAC') NCALCU = 1      
C
C=======================================================================
C    
C     INITIAL CONDITION FOR CONSTANT FLOW DISCHARGE
C
      IF(LCONDIS) THEN
        SISYPHE_CFD = LCONDIS
        NSIS_CFD    = NCONDIS
        CONSTFLOW   = .FALSE.
      ELSE
        SISYPHE_CFD = .FALSE.
        NSIS_CFD    = 1
        CONSTFLOW   = .FALSE.
      ENDIF
C
C=======================================================================
C
C     FIN DES INITIALISATIONS
         IF(DEBUG.GT.0) WRITE(LU,*) 'END_INITIALIZATION'
      ENDIF ! IF (PART==0 OR PART = -1)
C 
C=======================================================================
C
      IF(PART==1.OR.PART==-1) THEN
C 
        IF(DEBUG.GT.0) WRITE(LU,*) 'TIME_LOOP'
C
C=======================================================================
C
C : 3                    /* BOUCLE EN TEMPS */
C
C=======================================================================
C
C----   ARRET DU CALCUL SI LE NOMBRE D'ITERATIONS DEMANDE EST NUL
C
        IF(NIT == 0) THEN
          IF (LNG == 1) WRITE(LU,200)
          IF (LNG == 2) WRITE(LU,201)
200       FORMAT(' ARRET DANS SISYPHE, NOMBRE D''ITERATIONS',/,
     *         ' DEMANDE NUL')
201       FORMAT(' STOP IN SISYPHE, NUMBER OF ITERATIONS EQ.0')
          CALL PLANTE(1)
          STOP
        ENDIF    
C
C---------------------------------------------------------------------
C             DEBUT DES CALCULS
C---------------------------------------------------------------------
C       BOUCLE SUR LE NOMBRE D'EVENEMENTS 
C       (EN PERMANENT BOUCLE SUR LES PAS DE TEMPS)
C-----------------------------------
C
        IF(CODE(1:7) == 'TELEMAC') THEN
C         VALNIT IS THE NUMBER OF SISYPHE ITERATIONS UPDATING FOR COUPLING
          VALNIT = (TELNIT/PERICOU)*PERICOU-PERICOU+1
C         MODIF JMH+CV POUR EVITER D'APPELER VALIDA DEUX FOIS
C         QUAND ON APPELLE SUCCESSIVEMENT POUR CHARRIAGE ET SUSPENSION
          IF(GRAFCOUNT.GT.TELNIT) VALNIT=NIT+1         
C         CHARR, SUSP AND TIME STEP MONITORED BY CALLING PROGRAM
          CHARR = CHARR_TEL  
          SUSP= SUSP_TEL 
          AT0=T_TEL         
        ENDIF
C
        DO 710 MN = 1, NCALCU
C
          IF(.NOT.PERMA.AND.CODE(1:7).NE.'TELEMAC') THEN
C           DETERMINATION DU PREMIER ENREGISTREMENT A LIRE :
C           NUMDEB EST LE PREMIER ENREGISTREMENT A LIRE DANS LE FICHIER 
C           HYDRODYNAMIQUE.
            NUMDEB=NUMEN0
            IF(NUMDEB+NIDT > NUMEN) THEN
              IF (LNG == 1) WRITE(LU,202)
              IF (LNG == 2) WRITE(LU,203)
202           FORMAT(1X,'FICHIER HYDRODYNAMIQUE PAS ASSEZ LONG')
203           FORMAT(1X,'THE HYDRODYNAMIC FILE IS NOT LONG ENOUGH')             
              CALL PLANTE(1)
            ENDIF  
          ENDIF
C
C         BOUCLE SUR LES ENREGISTREMENTS  (SI PERMA NIDT=1)
C         ------------------------------
C         EN COUPLAGE ON NE FAIT QU'UN SEUL PAS DE LA BOUCLE
C  
CV          IF(CODE(1:7) == 'TELEMAC') DT=DT_TEL
C
          DO 700 MT = 1, NIDT
C
C  ----   DETERMINATION DU NUMERO DU PAS DE TEMPS :
C
            LT = (MN-1)*NIDT +  MT
C     
            IF(CODE(1:7) == 'TELEMAC') THEN 
              DT=DT_TEL
              LT    = LOOPCOUNT
              LEOPR = GRAFCOUNT
              LISPR = LISTCOUNT
              NSOUS=1
            ENDIF
C
C  ----    LOGIQUE DECIDANT DES IMPRESSIONS SUR LISTING :
C
            ENTETS = .FALSE.
            IF(LISPR*((LT-1+PERICOU)/LISPR) == (LT-1+PERICOU)) THEN
              ENTET = .TRUE.
            ELSE
              ENTET = .FALSE.
            ENDIF
C
C----      LECTURE ET REACTUALISATION DE H ET ZF
C----      SI 1ER PASSAGE OU NON PERMANENT ET NON COUPLAGE
C        
           DTS = DT / NSOUS
C
           ISOUS = 0
C
           IF(.NOT.PERMA.OR.PASS) THEN
C
C            ATTENTION : LA VALEUR DE ESOMT DU FICHIER SISHYD N'EST PAS LUE
C            NOTE : NOM DE SIS HYD SET TO ' ' IF COUPLING
C
             IF(SIS_FILES(SISHYD)%NAME(1:1).NE.' ')  THEN
C
C              WORK ON ZF,QU,QV,Z WILL BE IN FACT DONE ON:
C              T4,DEL_QU,DEL_QV AND DEL_Z
C              BY PLAYING WITH POINTERS              
               SAVEZF=>ZF%R
               SAVEQU=>QU%R
               SAVEQV=>QV%R
               SAVEZ =>Z%R
               ZF%R  =>T4%R
               QU%R  =>DEL_QU%R
               QV%R  =>DEL_QV%R
               Z%R   =>DEL_Z%R
C
               NUMDEB=NUMDEB+1
C
               IF(ENTET) WRITE(LU,*) 'DEFINITION INITIALE DES VITESSES'            
C
               CALL BIEF_SUITE(VARSOR,VARCL,NUMDEB,
     &            SIS_FILES(SISHYD)%LU,SIS_FILES(SISHYD)%FMT,
     &            HIST,0,NPOIN,BID,TEXTPR,VARCLA,0,
     &            TROUVE,ALIRE,ENTET,PERMA,MAXVAR)
C
               IF(DEBUG.GT.0) WRITE(LU,*) 'RESCUE_SISYPHE_NOTPERMA'
               CALL RESCUE_SISYPHE_NOTPERMA
     *              (QU%R,QV%R,Q%R,U2D%R,V2D%R,HN%R,Z%R,T4%R,
     *               HW%R,TW%R,THETAW%R,NPOIN,TROUVE,ALIRE,ICF,ENTET)
               IF(DEBUG.GT.0) WRITE(LU,*) 'END_RESCUE_SISYPHE_NOTPERMA'
C
C              BACK TO ORIGINAL ADDRESSES           
               ZF%R=>SAVEZF
               QU%R=>SAVEQU
               QV%R=>SAVEQV
               Z%R=>SAVEZ
C               
C              INCREMENT OF QU, QV AND Z PER SUB-TIME-STEP
               DO I = 1,NPOIN
                 DEL_QU%R(I) = (DEL_QU%R(I)-QU%R(I))/NSOUS
                 DEL_QV%R(I) = (DEL_QV%R(I)-QV%R(I))/NSOUS
                 DEL_Z%R(I)  = (DEL_Z%R(I) -Z%R(I)) /NSOUS
               ENDDO
C  
C              REACTUALISATION HYDRO NON PERMANENT
C              (A DEPLACER DANS RESCUE_SISYPHE_NOTPERMA)
C              -----------------------------------
C              CLIPPING OF NEGATIVE DEPTHS
C              CALCUL DE U2D ET V2D
C
               CALL OS('X=Y-Z   ', X=HN, Y=Z, Z=ZF)
C
               IF(OPTBAN.GT.0) THEN         
                 DO I = 1,NPOIN
                   IF(HN%R(I).LT.HMIN) THEN
                     U2D%R(I)=0.D0
                     V2D%R(I)=0.D0
                     HN%R(I) = MAX(HN%R(I),HMIN)     
                   ELSE
                     U2D%R(I)=QU%R(I)/HN%R(I)
                     V2D%R(I)=QV%R(I)/HN%R(I)
                   ENDIF
                 ENDDO
               ELSE
                 CALL OS('X=Y/Z   ', X=U2D, Y=QU,   Z=HN)
                 CALL OS('X=Y/Z   ', X=V2D, Y=QV,   Z=HN)                       
               ENDIF  
C  
               IF(DEBUG.GT.0) WRITE(LU,*) 'CONDIM_SISYPHE'
               CALL CONDIM_SISYPHE
     *        (U2D%R,V2D%R,QU%R,QV%R,HN%R,ZF%R,Z%R,ESOMT%R,THETAW%R,
     *         Q%R,HW%R,TW%R,MESH%X%R,MESH%Y%R,NPOIN,AT0,PMAREE)
               IF(DEBUG.GT.0) WRITE(LU,*) 'END_CONDIM_SISYPHE'
             ENDIF ! (SIS_FILES(SISHYD)%NAME(1:1) /=' ')
           ENDIF ! (NOT.PERMA.OR.PASS)
C                
        IF(PASS) THEN
C         LOGIQUE DE LECTURE MIS A FALSE EN PERMANENT
          IF (PERMA) PASS = .FALSE.
        ELSE
C         CALCUL DE LA HAUTEUR D'EAU
          CALL OS('X=Y-Z   ', X=HN, Y=Z, Z=ZF)
        ENDIF
C   
C       COUPLING
C
        IF(CODE(1:7) == 'TELEMAC') THEN
C
C         OV INSTEAD OF OS IN ORDER TO AVOID PROBLEMS WITH QUASI-BUBBLE ELEMENTS
C         ONE OPERATES ONLY ON THE (1:NPOIN) RANGE OF THE TELEMAC FIELDS
C         IT IS A HIDDEN DISCRETISATION CHANGE
C
          CALL OV( 'X=Y     ',U2D%R, U_TEL%R, U_TEL%R, 0.D0, NPOIN)
          CALL OV( 'X=Y     ',V2D%R, V_TEL%R, V_TEL%R, 0.D0, NPOIN)
          CALL OV( 'X=Y     ', HN%R, H_TEL%R, H_TEL%R, 0.D0, NPOIN)
C         ADDED BY JMH 01/07/2004 (ZF MAY BE MODIFIED BY CALLING PROGRAM)
          CALL OS('X=Y     ', X=ZF, Y=ZF_SIS)
C         CLIPPING OF NEGATIVE DEPTHS
          IF(OPTBAN.GT.0) THEN
            DO I = 1,HN%DIM1            
              IF(HN%R(I).LT.HMIN) THEN
               U2D%R(I)=0.D0
               V2D%R(I)=0.D0
               HN%R(I)=HMIN
              ENDIF
            ENDDO
          ENDIF
C         FREE SURFACE
          CALL OS('X=Y+Z   ', X=Z, Y=ZF, Z=HN)         
C         PRODUCT H*
!         CALL OS('X=YZ    ', X=QU, Y=U2D, Z=HN)
C         PRODUCT H*V
!         CALL OS('X=YZ    ', X=QV, Y=V2D, Z=HN)
C         DISCHARGE 
!         CALL OS('X=N(Y,Z)', X=Q, Y=QU, Z=QV)
C
        ENDIF
C
C       END OF COUPLING
! =========================================================================
! POUR TRAITEMENT BANCS DECOUVRANTS, DEFINITION DES MASQUES 
! =====================================================================!
C               
        IF(OPTBAN.EQ.2) THEN
C
C ----    CREATION D'UN TABLEAU DE MASQUAGE PAR ELEMENTS
C
          CHGMSK = .FALSE.
          CALL OS ('X=Y     ', X=MSKTMP, Y=MASKEL)
          CALL OS ('X=C     ', X=MASKEL, C=1.D0)
          IF(CODE(1:7) == 'TELEMAC') THEN
C           ON FABRIQUE LES MASQUES A PARTIR DES VALEURS DE H
C           NON CLIPPEES FOURNIES PAR TELEMAC
            CALL MASKTF(MASKEL%R,H_TEL%R,HMIN,MESH%IKLE%I,
     *                  NELEM,NPOIN)
          ELSE
            CALL MASKTF(MASKEL%R,HN%R,HMIN,MESH%IKLE%I,
     *                  NELEM,NPOIN)
          ENDIF  
C
          DO I=1,NELEM
            IF(MASKEL%R(I).NE.MSKTMP%R(I)) THEN
              CHGMSK = .TRUE.
              EXIT 
            ENDIF
          ENDDO
C
C        JMH 17/12/2009
C
C        ELSEIF(OPTBAN.EQ.1) THEN
C
C          CANCELS Q QU AND QV IF HN.LE.0.D0
C          CALL MASKAB(HN%R,Q%R,QU%R,QV%R,NPOIN)
C     
        ENDIF
C
C ----   CONSTRUCTION DU MASQUE DES POINTS A PARTIR DU MASQUE DES ELEMENTS
C ----   ET CHANGEMENT DE IFAMAS (IFABOR AVEC MASQUAGE)
C
        IF(MSK) CALL MASKTO(MASKEL%R,MASKPT,IFAMAS%I,
     *                      MESH%IKLE%I,
     *                      MESH%IFABOR%I,MESH%ELTSEG%I,MESH%NSEG,
     *                      NELEM,NPOIN,IELMT,MESH)
C
C ------------------------------------------------------------------
C  DEBUT SOUS-ITERATIONS EN NON-PERMANENT
C
C ------------------------------------------------------------------
C 
702      CONTINUE
C
         ISOUS = ISOUS + 1
         AT0=AT0+DTS
         IF(ENTET.AND.ISOUS.EQ.1) CALL ENTETE_SISYPHE(2,AT0,LT)
         IF(ENTET.AND.ISOUS.EQ.NSOUS) ENTETS=.TRUE.
C                   
C---------------------------------------------------------------------
C              COEFFICIENT DE FROTTEMENT VARIABLE EN TEMPS
C---------------------------------------------------------------------
C          
         CALL CORSTR_SISYPHE
C
C ----   LECTURE DES CONDITIONS AUX LIMITES 
C
         CALL CONLIT(MESH%NBOR%I) 
C
C =======================================================================
C
C        SI "PAS DE TEMPS VARIABLE" = OUI NSOUS SERA DETERMINE PLUS LOIN
C        LE CALCUL DU PAS DE TEMPS SIS EST DEPALCE AVANT LA LECTURE DES 
C        CONDITIONS HYDRO
C
C  ---   DIAMETRE MOYEN DE LA COUCHE ACTIVE ET DE LA SOUS-COUCHE
C  
         IF(.NOT.MIXTE.AND.NSICLA.GT.1) CALL MEAN_GRAIN_SIZE
C         
C --- CALCUL DE LA VITESSE MOYESSE UNORM 
C
         CALL OS('X=N(Y,Z)',X=UNORM,Y=U2D,Z=V2D)
C         
C  ---     CALCUL DE LA VITESSE ORBITALE DE LA HOULE --> UW
C
         IF(HOULE) THEN
           CALL CALCUW(UW%R,HN%R,HW%R,TW%R,GRAV,NPOIN)
         ENDIF
C
         CALL TOB_SISYPHE
     &     (TOB,TOBW, MU, KS, KSP,KSR,CF, FW, 
     &      CHESTR, UETCAR, CF_TEL, CODE ,
     &      KFROT, ICR, KSPRATIO,HOULE,
     &      GRAV,XMVE,  XMVS, VCE, KARMAN,ZERO,
     &      HMIN,HN, ACLADM, UNORM,UW, TW, NPOIN)
C
C  FIN MODULE D'INITIALISASTION
C 
      ! ******************** !
      ! BED-LOAD COMPUTATION !
      ! ******************** !

        IF(CHARR) THEN 
          IF(DEBUG.GT.0) WRITE(LU,*) 'BEDLOAD_MAIN'
          CALL BEDLOAD_MAIN
     &        (ACLADM,KSP,KSR,VOLU2D,UNSV2D, 
     &         CF,EBOR,FW,HN,LIQBOR,MASK,MASKEL,
     &         MASKPT,Q,QBOR,U2D,V2D,S,UNLADM,UW,THETAW,
     &         MU,TOB,TOBW,TW,ZF,DEBUG,HIDFAC,ICF,
     &         IELMT,ISOUS,KDDL,KDIR,KENT,KINC,KLOG,KNEU,KSORT,
     &         LOADMETH,LT,NPOIN,NPTFR,NSICLA,
     &         OPTBAN,LS0,BETA,FD90,FDM,GRAV,HIDI,HMIN,
     &         VCE,XKV,XMVE,XMVS,XWC,PI,KARMAN,ZERO,
     &         KARIM_HOLLY_YANG,MSK,SUSP,VF,ENTET,
     &         CONST_ALAYER,LCONDIS,LGRAFED,MESH,
     &         ELAY,LIEBOR,LIMTEC,MASKTR,
     &         IT1,T1,T2,T3,T4,T5,T6,T7,T8,T9,
     &         T10,T11,T12,T13,UNORM,AC,AT0,DTS,ELAY0,FRACSED_GF,
     &         AVAIL,BREACH,CALFA,COEFPN,DZF_GF,HIDING,
     &         QSCL_C,QSCL_S,QS_C,QSCLXC,QSXC,QSCLYC,
     &         QSYC,SALFA,ZF_C,ZFCL_C,NSOUS,ENTETS,
     &         SECCURRENT,SLOPEFF,PHISED,DEVIA,BETA2,BIJK,
     &         SEDCO,HOULE)
C
          IF(DEBUG.GT.0) WRITE(LU,*) 'END_BEDLOAD_MAIN'
C 
C         UPDATING BOTTOM
C
          CALL OS('X=X+Y   ',X=ZF,Y=ZF_C) 
C
C         UPDATING LAYERS  --> ELAY
C
          IF(.NOT.MIXTE.AND.NSICLA.GT.1) THEN
            IF(DEBUG.GT.0) WRITE(LU,*) 'LAYER'
            CALL LAYER(ZFCL_C,NLAYER,ZR,ZF,ESTRAT,ELAY,VOLU2D,
     &                 ACLADM,NSICLA,NPOIN,ELAY0,VOLTOT,ES,
     &                 AVAIL,CONST_ALAYER,DTS,T2%R,IT1%I)
            IF(DEBUG.GT.0) WRITE(LU,*) 'END_LAYER'
          ELSE
            CALL OS('X=Y-Z   ',X=ELAY,Y=ZF,Z=ZR)
          ENDIF
C fin charriage
        ENDIF

      ! ********************** !
      ! SUSPENSION COMPUTATION ! 
      ! ********************** !
        IF(SUSP) THEN
C           
          IF(DEBUG.GT.0) WRITE(LU,*) 'SUSPENSION_MAIN'
          CALL SUSPENSION_MAIN
     &(SLVTRA,HN,HN_TEL,MU,TOB,ACLADM,KSP,KSR,CF,
     * VOLU2D,V2DPAR,UNSV2D,AFBOR,BFBOR,ZF,LICBOR,
     * IFAMAS,MASKEL,MASKPT,U2D,V2D,NSICLA,
     & NPOIN,NPTFR,IELMT,OPTDIF,RESOL,LT,NIT,OPTBAN,OPTSUP,
     & OPDTRA,KENT,KSORT,KLOG,KINC,KNEU,KDIR,KDDL,ISOUS,NSOUS,
     & DEBUG,DTS,CSF_VASE,CSF_SABLE,ZERO,GRAV,XKX,XKY,
     & KARMAN,XMVE,XMVS,HMIN,XWC,VITCD,VITCE,PARTHENIADES, 
     & ENTET,BILMA,MSK,CHARR,IMP_INFLOW_C,MESH,ZF_S,CS,
     & CST,CTILD,CBOR,DISP,IT1,IT2,IT3,IT4,TB,T1,T2,T3,T4,T5,T6,
     & T7,T8,T9,T10,T11,T12,W1,TE1,CLT,TE2,TE3,S,AM1_S,AM2_S,MBOR,
     & ELAY,LIMDIF,MASKTR,TETA_SUSP,AC,
     & MASED0,MASINI,MASTEN,MASTOU,ES,AVAIL,ENTETS,PASS_SUSP,
     & ZFCL_S,HPROP,FLUDPT,FLUDP,FLUER,DISP_C,KX,KY,KZ,UCONV,  
     & VCONV,QSXS,QSYS,QSCLXS,QSCLYS,QSCL_S,QS_S,QS_C,
     & CSTAEQ,ICQ,MASTCP,MASFIN,MASDEPT,MASDEP,MASSOU,CORR_CONV,
     & ZREF,SEDCO,VISC_TEL,CODE,DIFT,DM1,UCONV_TEL,VCONV_TEL,
     & ZCONV,SOLSYS,FLBOR_TEL,FLBOR_SIS,FLBORTRA,NUMLIQ%I,NFRLIQ,
     & MIXTE,NCOUCH_TASS,CONC_VASE,TOCE_VASE,
     & FLUER_VASE,TOCE_MIXTE,MS_SABLE,MS_VASE,TASS) 
         IF(DEBUG.GT.0) WRITE(LU,*) 'END_SUSPENSION_MAIN'
C
C      UPDATING BOTTOM
C
       CALL OS('X=X+Y   ',X=ZF,Y=ZF_S) 
C
C      UPDATING LAYERS
C      Redinition de la couche de sediment erodable
C      granulo etendue (a remplacer par NOMBLAY>1
C
        IF(.NOT.MIXTE.AND.NSICLA.GT.1) THEN
          IF(DEBUG.GT.0) WRITE(LU,*) 'LAYER'
          CALL LAYER(ZFCL_S,NLAYER,ZR,ZF,ESTRAT,ELAY,VOLU2D,
     &               ACLADM,NSICLA,NPOIN,ELAY0,VOLTOT,ES,
     &               AVAIL,CONST_ALAYER,DTS,T2%R,IT1%I)
          IF(DEBUG.GT.0) WRITE(LU,*) 'END_LAYER'
        ELSE  
          CALL OS('X=Y-Z   ',X=ELAY,Y=ZF,Z=ZR)
        ENDIF    
C fin suspension
      ENDIF
!
! RECONSTITUTION DES DONNEES CHARRIAGE ET/OU SUSPENSION
! -----------------------------------------------------
!
        IF( DEBUG.GT.0) WRITE(LU,*) 'QS_RESULT'
!        
        CALL OS('X=0     ', X=QSX)
        CALL OS('X=0     ', X=QSY)
!
        DO I = 1, NSICLA
          CALL OS('X=Y+Z   ', X=QSCLX%ADR(I)%P, Y=QSCLXC%ADR(I)%P,
     &                                          Z=QSCLXS%ADR(I)%P)
          CALL OS('X=Y+Z   ', X=QSCLY%ADR(I)%P, Y=QSCLYC%ADR(I)%P,
     &                                          Z=QSCLYS%ADR(I)%P)
          CALL OS('X=N(Y,Z)', X=QSCL%ADR(I)%P,  Y=QSCLX%ADR(I)%P,
     &                                          Z=QSCLY%ADR(I)%P)
          CALL OS('X=X+Y   ', X=QSX, Y=QSCLX%ADR(I)%P)
          CALL OS('X=X+Y   ', X=QSY, Y=QSCLY%ADR(I)%P)
        ENDDO
        CALL OS('X=N(Y,Z)', X=QS, Y=QSX, Z=QSY)
        IF(DEBUG.GT.0) WRITE(LU,*) 'END_QS_RESULT'
!
!=======================================================================
!
!     PENTE MAXIMUM DU FOND: evol mis dans T1
!
      IF(SLIDE) THEN
!
        IF(ENTET) CALL ENTETE_SISYPHE(14,AT0,LT)
!
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE MAXSLOPE'
        CALL MAXSLOPE(PHISED,ZF%R,ZR%R,MESH%XEL%R,MESH%YEL%R,MESH%NELEM,
     *                MESH%NELMAX,NPOIN,MESH%IKLE%I,T1,UNSV2D,MESH)
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE MAXSLOPE'
        CALL OS('X=X+Y   ',X=ZF,Y=T1) 
! 
        IF(NSICLA.EQ.1) THEN
          CALL OS('X=Y-Z   ',X=ELAY,Y=ZF,Z=ZR)
        ELSE
          WRITE(LU,*) 'SLIDE NOT IMPLEMENTED WITH GRADED SEDIMENT'
          CALL PLANTE(1)
          STOP
        ENDIF       
!
      ENDIF
!
!========================================================================
!
!     SETTLING: EVOLUTION COMPUTED IN T3
!
      IF(TASS) THEN
C
        IF(ENTET) THEN
          IF(.NOT.CHARR.AND..NOT.SUSP.AND..NOT.SLIDE) THEN
            CALL ENTETE_SISYPHE(2,AT0,LT)
          ENDIF
          CALL ENTETE_SISYPHE(15,AT0,LT)
        ENDIF
C
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE TASSEMENT'
        CALL TASSEMENT(ZF,NPOIN,DTS,ELAY,T3,T2,LT,AVAIL,NSICLA,
     *                 ES,XMVS,XKV,TRANS_MASS,CONC_VASE,NCOUCH_TASS, 
     *                 MS_SABLE%R,MS_VASE%R)        
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE TASSEMENT'
C
C       UPDATING ZF (ELAY HAS BEEN UPDATED IN TASSEMENT)
C
        CALL OS('X=X+Y   ',X=ZF,Y=T3) 
C        
      ENDIF     
!
!=======================================================================
! : 5        CALCUL DES EVOLUTIONS SUR CE CYCLE DE PAS DE TEMPS
!            ET REACTUALISATION APRES CE CALCUL
!=======================================================================
!
! ----  CALCUL DES EVOLUTIONS SUR CE (SOUS) PAS DE TEMPS
!
      IF(CHARR) THEN
        CALL OS('X=Y     ',X=E,Y=ZF_C)
      ELSE
        CALL OS('X=0     ',X=E)
      ENDIF
      IF(SUSP)  CALL OS('X=X+Y   ',X=E,Y=ZF_S)
      IF(SLIDE) CALL OS('X=X+Y   ',X=E,Y=T1) 
      IF(TASS)  CALL OS('X=X+Y   ',X=E,Y=T3)
!
      CALL OS('X=X+Y   ', X=ESOMT, Y=E)
C
C  REACTUALISATIONS 
C    
        IF(PART.EQ.-1) THEN
C
        CALL OS('X=X-Y   ',X=HN,Y=E)
        IF(OPTBAN.GT.0) THEN
          DO I = 1,HN%DIM1
            IF(HN%R(I).LT.HMIN) THEN
              U2D%R(I)=0.D0
              V2D%R(I)=0.D0
              HN%R(I) =HMIN
            ELSE
              U2D%R(I)= QU%R(I)/HN%R(I)
              V2D%R(I)= QV%R(I)/HN%R(I)
            ENDIF  
          ENDDO
        ELSE
          CALL OS('X=Y/Z   ', X=U2D, Y=QU,   Z=HN)
          CALL OS('X=Y/Z   ', X=V2D, Y=QV,   Z=HN)                       
        ENDIF
C
C=======================================================================
C : 6     ARRET SI EVOLUTIONS SUPERIEURES A EMAX = RC*(HAUTEUR INITIALE)
C=======================================================================
C
C       DETERMINATION DU SEUIL D'EVOLUTION MAXIMUM
        DO I = 1, NPOIN
          EMAX%R(I) = RC*MAX(HN%R(I),HMIN)
        ENDDO
C 
C ----  TEST D'ARRET LORSQUE LES EVOLUTIONS DEPASSENT UN CERTAIN SEUIL
C CCV CE TEST N'INTERVIENT QUE EN MODE SISYPHE SEUL
C
        IF(DEBUG.GT.0) WRITE(LU,*) 'ARRET'
          CALL SIS_ARRET(ESOMT,EMAX,HN,VARSOR,NPOIN,MN,
     *                   SIS_FILES(SISRES)%LU,SIS_FILES(SISRES)%FMT,
     *                   MAXVAR,AT0,RC,HIST,BINRESSIS,TEXTE,
     *                   SORLEO,SORIMP,T1,T2)
          IF(DEBUG.GT.0) WRITE(LU,*) 'END_ARRET'
        ENDIF
C
C ----     CONSTANT FLOW DISCHARGE
C     
        IF(LCONDIS) THEN
          CALL CONDIS_SISYPHE(CONSTFLOW)
        ELSE
          CONSTFLOW =.FALSE.
        ENDIF
C
C=======================================================================
C : 8     BILAN DE MASSE
C=======================================================================
C       CALCUL DES COMPOSANTES DU DEBIT SOLIDE POUR LE BILAN DE MASSE, LES
C       SORTIES GRAPHIQUES ET L'ETAPE DE VALIDATION
C
        IF(BILMA.AND.CHARR) THEN
          IF(DEBUG.GT.0) WRITE(LU,*) 'BILAN_SISYPHE'
          CALL BILAN_SISYPHE(E,ESOMT,QSX,QSY,
     &                       MESH,MSK,MASKEL,T1,T2,S,IELMU_SIS,VCUMU,
     &       DTS,NPTFR,ENTETS,ZFCL_C,QSCLXC,QSCLYC,NSICLA, 
     &       VOLTOT,DZF_GF,MASS_GF,LGRAFED,NUMLIQ%I,NFRLIQ)
          IF(DEBUG.GT.0) WRITE(LU,*) 'END_BILAN_SISYPHE'
        ENDIF
C
C       SECTIONS DE CONTROLE
C
        IF(NCP.GT.0) THEN
          IF(DEBUG.GT.0) WRITE(LU,*) 'FLUSEC_SISYPHE'
          CALL FLUSEC_SISYPHE(U2D,V2D,HN,
     *                        QSXC,QSYC,CHARR,QSXS,QSYS,SUSP,
     *                        MESH%IKLE%I,
     *                        MESH%NELMAX,MESH%NELEM,
     *                        MESH%X%R,MESH%Y%R,
     *                        DT,NCP,CTRLSC,ENTETS,AT0,MESH%KNOGL%I)
          IF(DEBUG.GT.0) WRITE(LU,*) 'END_FLUSEC_SISYPHE'
        ENDIF      
C
C-----------------------------------------------------------------------
C
        IF(.NOT.PERMA.AND.SIS_FILES(SISHYD)%NAME(1:1).NE.' ') THEN
C 
C         REACTUALISATION HYDRODYNAMIQUE
C
C         IF READING ON HYDRODYNAMIC FILE, INCREMENTING QU, QV AND Z
          IF(SIS_FILES(SISHYD)%NAME(1:1).NE.' ')  THEN
            CALL OS('X=X+Y   ', X=QU, Y=DEL_QU)
            CALL OS('X=X+Y   ', X=QV, Y=DEL_QV)
            CALL OS('X=X+Y   ', X=Z , Y=DEL_Z)
          ENDIF
          CALL OS('X=Y-Z   ', X=HN, Y=Z, Z=ZF)
C         CLIPPING OF NEGATIVE DEPTHS
          IF(OPTBAN.GT.0) THEN
            DO I = 1, NPOIN
             IF(HN%R(I).LT.HMIN) THEN
               U2D%R(I)=0.D0
               V2D%R(I)=0.D0
               HN%R(I) = MAX(HN%R(I),HMIN)
             ELSE
               U2D%R(I)= QU%R(I)/HN%R(I)
               V2D%R(I)= QV%R(I)/HN%R(I)
             ENDIF  
            ENDDO
          ELSE
            CALL OS('X=Y/Z   ', X=U2D, Y=QU,   Z=HN)
            CALL OS('X=Y/Z   ', X=V2D, Y=QV,   Z=HN)                       
          ENDIF
C
        ENDIF
C                  
C       FIN DE LA BOUCLE SUR LES NSOUS SOUS PAS DE TEMPS
C ---------------------------------------------------------
        IF(DEBUG.GT.0) WRITE(LU,*) 'SOUS_ITERATION_NEXT'
        IF (ISOUS < NSOUS) GOTO 702
        IF(DEBUG.GT.0) WRITE(LU,*) 'END_SOUS_ITERATION'
C=======================================================================
C : 9        IMPRESSIONS DES VALEURS EXTREMES
C=======================================================================
C
        IF(ENTET.AND.CHARR) THEN
          WRITE(LU,*)
          CALL MAXI(XMAX,IMAX,E%R,NPOIN)
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
          IF(LNG.EQ.1) WRITE(LU,371) XMA,IMA
          IF(LNG.EQ.2) WRITE(LU,372) XMA,IMA
371       FORMAT(' EVOLUTION MAXIMUM        : ',G16.7,' NOEUD : ',I6)
372       FORMAT(' MAXIMAL EVOLUTION        : ',G16.7,' NODE  : ',I6)
          CALL MINI(XMIN,IMIN,E%R,NPOIN)
          IF(NCSIZE.GT.1) THEN
            XMI=P_DMIN(XMIN)
            IF(XMIN.EQ.XMI) THEN
              IMI=MESH%KNOLG%I(IMIN)
            ELSE
              IMI=0
            ENDIF
            IMI=P_IMAX(IMI)
          ELSE
            IMI=IMIN
            XMI=XMIN 
          ENDIF
          IF(LNG.EQ.1) WRITE(LU,373) XMI,IMI
          IF(LNG.EQ.2) WRITE(LU,374) XMI,IMI
373       FORMAT(' EVOLUTION MINIMUM        : ',G16.7,' NOEUD : ',I6)
374       FORMAT(' MINIMAL EVOLUTION        : ',G16.7,' NODE  : ',I6)
C
          IF(CONST_ALAYER) THEN
            IF(NSICLA.GT.1.AND.XMI.LT.-0.5D0*ELAY0) THEN
              IF(LNG.EQ.1) WRITE(LU,885) 
              IF(LNG.EQ.2) WRITE(LU,886) 
885           FORMAT(' EROSION SUPERIEURE A EPAISSEUR DE COUCHE !')
886           FORMAT(' EROSION GREATER THAN ONE LAYER THICKNESS !')
            ENDIF
            IF(NSICLA.GT.1.AND.XMA.GT.ELAY0) THEN
              IF(LNG.EQ.1) WRITE(LU,887) 
              IF(LNG.EQ.2) WRITE(LU,888) 
887           FORMAT(' DEPOT SUPERIEUR A EPAISSEUR DE COUCHE !')
888           FORMAT(' DEPOSITION MORE THAN ONE LAYER THICKNESS !') 
            ENDIF
          ELSE
            DO J=1,NPOIN
              ELAY0 = 3.D0*ACLADM%R(J)
              IF(NSICLA.GT.1.AND.E%R(J).LT.-0.5D0*ELAY0) THEN
                IF(LNG.EQ.1) WRITE(LU,885) 
                IF(LNG.EQ.2) WRITE(LU,886) 
              ENDIF
              IF(NSICLA.GT.1.AND.E%R(J).GT.ELAY0) THEN
                IF(LNG.EQ.1) WRITE(LU,887) 
                IF(LNG.EQ.2) WRITE(LU,888) 
              ENDIF
            ENDDO
          ENDIF
        ENDIF
        IF(ENTET) THEN
          CALL MAXI(XMAX,IMAX,ESOMT%R,NPOIN)
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
          IF (LNG.EQ.1) WRITE(LU,881) XMA,IMA
          IF (LNG.EQ.2) WRITE(LU,882) XMA,IMA
881       FORMAT(' EVOLUTION MAXIMUM TOTALE : ',G16.7,' NOEUD : ',I6)
882       FORMAT(' TOTAL MAXIMAL EVOLUTION  : ',G16.7,' NODE  : ',I6)
          CALL MINI(XMIN,IMIN,ESOMT%R,NPOIN)
          IF(NCSIZE.GT.1) THEN
            XMI=P_DMIN(XMIN)
            IF(XMIN.EQ.XMI) THEN
              IMI=MESH%KNOLG%I(IMIN)
            ELSE
              IMI=0
            ENDIF
            IMI=P_IMAX(IMI)
          ELSE
            IMI=IMIN
            XMI=XMIN 
          ENDIF
          IF (LNG.EQ.1) WRITE(LU,883) XMI,IMI
          IF (LNG.EQ.2) WRITE(LU,884) XMI,IMI
883       FORMAT(' EVOLUTION MINIMUM TOTALE : ',G16.7,' NOEUD : ',I6)
884       FORMAT(' TOTAL MINIMAL EVOLUTION  : ',G16.7,' NODE  : ',I6)
        ENDIF

C=======================================================================
C : 10         IMPRESSION DES RESULTATS A CE PAS DE TEMPS
C              ET COMPARAISON AVEC UN FICHIER DE REFERENCE
C=======================================================================
C       

C       IN COUPLING MODE, OUTPUT TIMES OF TELEMAC AND SISYPHE ARE
C       SYNCHRONISED, IT MEANS THAT WE MUST HAVE :
C       LT * DT (TIME AT THE END OF TIME STEP LT IN TELEMAC)
C       EQUAL TO:
C       (LT-1)*DT + PERICOU*DT (TIME AT THE END OF TIME STEP LT IN SISYPHE
C       HENCE THE FACT THAT LT IS REPLACED BY LT-1+PERICOU
C       DEFAULT VALUE OF PERICOU IS 1.
C
        IF(UNIT) CALL OS('X=CX    ', X=CS, C= XMVS)      
C
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE PREDES'
        CALL PREDES(LT-1+PERICOU,AT0)
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE PREDES'
        IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE BIEF_DESIMP'
        CALL BIEF_DESIMP(SIS_FILES(SISRES)%FMT,VARSOR,
     *                   HIST,0,NPOIN,SIS_FILES(SISRES)%LU,
     *                   'STD',AT0,LT-1+PERICOU,LISPR,LEOPR,
     *                   SORLEO,SORIMP,MAXVAR,TEXTE,PTINIG,PTINIL)
        IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE BIEF_DESIMP'
C
        IF(UNIT) CALL OS('X=CX    ', X=CS,  C= 1.D0/XMVS)      
C
C       SEND THE NEW ZF TO TELEMAC-2D OR 3D
C 
        IF(CODE(1:7) == 'TELEMAC') THEN     
          CALL OV ('X=Y     ', ZF_SIS%R, ZF%R, ZF%R, 0.D0, NPOIN)
        ENDIF
C
C       LE SOUS-PROGRAMME VALIDA DE LA BIBLIOTHEQUE EST STANDARD
C       IL PEUT ETRE MODIFIE POUR CHAQUE CAS PARTICULIER
C       L'APPEL DOIT ETRE LAISSE DANS LA BOUCLE EN TEMPS
C
        IF(VALID) THEN
          IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE BIEF_VALIDA'
          CALL BIEF_VALIDA(TB,TEXTPR,SIS_FILES(SISREF)%LU,
     *                     SIS_FILES(SISREF)%FMT,
     *                     VARSOR,TEXTE,SIS_FILES(SISRES)%LU,
     *                     SIS_FILES(SISRES)%FMT,
     *                     MAXVAR,NPOIN,LT,VALNIT,ALIRV)
          IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE BIEF_VALIDA'
        ENDIF
C
C       FIN DE LA BOUCLE SUR LES ENREGISTREMENTS : 700
700     CONTINUE
C
C=======================================================================
C
C       FIN DE LA BOUCLE SUR LE NOMBRE D'EVENEMENTS : 710
C
710     CONTINUE
C
C-----------------------------------------------------------------------
C
        IF(DEBUG.GT.0) WRITE(LU,*) 'END_TIME_LOOP'
      ENDIF
C
C-----------------------------------------------------------------------
C
      IF(DREDGESIM.AND.(LOOPCOUNT.EQ.TELNIT.AND.PART.EQ.1.
     &                                            .OR. PART.EQ.-1)) THEN
         CALL DREDGESIM_INTERFACE(3)
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
