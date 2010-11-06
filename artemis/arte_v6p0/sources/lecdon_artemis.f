C                       *************************
                        SUBROUTINE LECDON_ARTEMIS
C                       *************************
C
     * (FILE_DESC,PATH,NCAR,CODE)
C
C***********************************************************************
C
C  ARTEMIS VERSION 6.0      20/04/99  D. AELBRECHT (LNH) 01 30 87 74 12
C
C  LINKED TO BIEF VERS. 5.0  17/08/94  J-M HERVOUET (LNH) 01 30 87 80 18
C
C***********************************************************************
C
C  FUNCTION : READING THE PARAMETER FILE THROUGH A DAMOCLES CALL
C
C----------------------------------------------------------------------
C
C CALLED BY : HOMERE_ARTEMIS
C
C SUBROUTINE CALLED : DAMOC
C
C**********************************************************************
C
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_ARTEMIS   
      USE INTERFACE_ARTEMIS, EX_LECDON_ARTEMIS => LECDON_ARTEMIS 
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      LOGICAL LSTOP
      INTEGER I,KK
C
      CHARACTER*8 MNEMO(MAXVAR)
C
C-----------------------------------------------------------------------
C
C DECLARATION DES TABLEAUX POUR L'APPEL DE DAMOCLES
C
      INTEGER, PARAMETER :: NMAX = 300
C
      INTEGER ADRESS(4,NMAX),DIMENS(4,NMAX)
      DOUBLE PRECISION MOTREA(NMAX)
      INTEGER          MOTINT(NMAX)
      LOGICAL          MOTLOG(NMAX)
      CHARACTER*144    MOTCAR(NMAX)
      CHARACTER*72     MOTCLE(4,NMAX,2)
      INTEGER          TROUVE(4,NMAX) 
      LOGICAL DOC
      CHARACTER(LEN=250) :: NOM_CAS
      CHARACTER(LEN=250) :: NOM_DIC
!ARGUMENTS
      CHARACTER(LEN=24), INTENT(IN)     :: CODE
      CHARACTER(LEN=144), INTENT(INOUT) :: FILE_DESC(4,NMAX)
      INTEGER, INTENT(IN)               :: NCAR
      CHARACTER(LEN=250), INTENT(IN)    :: PATH
C            
C FIN DES VARIABLES POUR DAMOCLES :
C
C-----------------------------------------------------------------------
C
C INITIALISATIONS POUR APPEL A DAMOCLES :
C
      DO 10 KK=1,NMAX
C
C       UN FICHIER NON DONNE PAR DAMOCLES SERA RECONNU PAR UN BLANC
C       (IL N'EST PAS SUR QUE TOUS LES COMPILATEURS INITIALISENT AINSI)
C
        MOTCAR(KK)(1:1)=' '
C
        DIMENS(1,KK) = 0
        DIMENS(2,KK) = 0
        DIMENS(3,KK) = 0
        DIMENS(4,KK) = 0
C
10    CONTINUE
C
C     IMPRESSION DE LA DOC
C
      DOC = .FALSE.
C
C-----------------------------------------------------------------------
C     OUVERTURE DES FICHIERS DICTIONNAIRE ET CAS
C-----------------------------------------------------------------------
C
      IF(NCAR.GT.0) THEN
C
        NOM_DIC=PATH(1:NCAR)//'ARTDICO'
        NOM_CAS=PATH(1:NCAR)//'ARTCAS'
C
      ELSE
C
        NOM_DIC='ARTDICO'
        NOM_CAS='ARTCAS'
C
      ENDIF
C
      OPEN(2,FILE=NOM_DIC,FORM='FORMATTED',ACTION='READ')
      OPEN(3,FILE=NOM_CAS,FORM='FORMATTED',ACTION='READ')
C
      CALL DAMOCLE( ADRESS , DIMENS , NMAX   , DOC     , LNG    , LU ,
     *              MOTINT , MOTREA , MOTLOG , MOTCAR  , MOTCLE ,
     *              TROUVE , 2   , 3   , .FALSE. , FILE_DESC)
C
C     DECRYPTAGE DES CHAINES SUBMIT
C
      CALL READ_SUBMIT(ART_FILES,44,CODE,FILE_DESC,300)
C
C-----------------------------------------------------------------------
C
      DO I=1,44
        IF(ART_FILES(I)%TELNAME.EQ.'ARTGEO') THEN
          ARTGEO=I
        ELSEIF(ART_FILES(I)%TELNAME.EQ.'ARTCAS') THEN
          ARTCAS=I
        ELSEIF(ART_FILES(I)%TELNAME.EQ.'ARTCLI') THEN
          ARTCLI=I
        ELSEIF(ART_FILES(I)%TELNAME.EQ.'ARTFON') THEN
          ARTFON=I
        ELSEIF(ART_FILES(I)%TELNAME.EQ.'ARTRES') THEN
          ARTRES=I
        ELSEIF(ART_FILES(I)%TELNAME.EQ.'ARTREF') THEN
          ARTREF=I
        ELSEIF(ART_FILES(I)%TELNAME.EQ.'ARTBI1') THEN
          ARTBI1=I
        ELSEIF(ART_FILES(I)%TELNAME.EQ.'ARTBI2') THEN
          ARTBI2=I
        ELSEIF(ART_FILES(I)%TELNAME.EQ.'ARTFO1') THEN
          ARTFO1=I
        ELSEIF(ART_FILES(I)%TELNAME.EQ.'ARTFO2') THEN
          ARTFO2=I
        ELSEIF(ART_FILES(I)%TELNAME.EQ.'ARTRBI') THEN
          ARTRBI=I
        ELSEIF(ART_FILES(I)%TELNAME.EQ.'ARTRFO') THEN
          ARTRFO=I
        ENDIF
      ENDDO
C
C-----------------------------------------------------------------------
C
C AFFECTATION DES PARAMETRES SOUS LEUR NOM EN FORTRAN :
C
C-----------------------------------------------------------------------
C
C MOTS-CLES ENTIERS :
C
         LEOPRD    = MOTINT( ADRESS(1, 1) )
         LISPRD    = MOTINT( ADRESS(1, 2) )
         SLVART%NITMAX    = MOTINT( ADRESS(1,3) )
         SLVART%PRECON    = MOTINT( ADRESS(1,4) )
         DISESP    = MOTINT( ADRESS(1, 5) )
         STDGEO    = MOTINT( ADRESS(1, 6) )
         STDRES    = MOTINT( ADRESS(1, 7) )
         SLVART%SLV = MOTINT( ADRESS(1,8) )
         LISFON    = MOTINT( ADRESS(1, 9) )
         NPALE     = MOTINT( ADRESS(1,10) )
         NDALE     = MOTINT( ADRESS(1,11) )
         IBREAK    = MOTINT( ADRESS(1,12) )
         NITDIS    = MOTINT( ADRESS(1,13) )
         REGIDO    = MOTINT( ADRESS(1,14) )
         FORMFR    = MOTINT( ADRESS(1,15) )
         SLVART%KRYLOV = MOTINT( ADRESS(1,16) )
         LVMAC     = MOTINT( ADRESS(1,17) )
         KFROT     = MOTINT( ADRESS(1,18) )
         OPTASS    = MOTINT( ADRESS(1,19) )
         PRODUC    = MOTINT( ADRESS(1,20) )
         MARDAT(1) = MOTINT( ADRESS(1,21) )
         MARDAT(2) = MOTINT( ADRESS(1,21) + 1 )
         MARDAT(3) = MOTINT( ADRESS(1,21) + 2 )
         MARTIM(1) = MOTINT( ADRESS(1,22) )
         MARTIM(2) = MOTINT( ADRESS(1,22) + 1 )
         MARTIM(3) = MOTINT( ADRESS(1,22) + 2 )
         NPRIV     = MOTINT( ADRESS(1,23) )
         NCSIZE    = MOTINT( ADRESS(1,24) )
C        ORIGIN COORDINATES
         I_ORIG    = MOTINT( ADRESS(1,25)   )
         J_ORIG    = MOTINT( ADRESS(1,25)+1 )
C POUR L'INSTANT ON DEBUTE LES IMPRESSIONS A ZERO
         PTINIG    = 0
         PTINIL    = 0
C
C MOTS-CLES REELS :
C
         PER       = MOTREA( ADRESS(2, 1) )
         TETAH     = MOTREA( ADRESS(2, 2) )
         GRAV      = MOTREA( ADRESS(2, 3) )
         SLVART%ZERO = MOTREA( ADRESS(2, 4) )
         SLVART%EPS = MOTREA( ADRESS(2, 5) )
         HMIN      = MOTREA( ADRESS(2, 6) )
         COTINI    = MOTREA( ADRESS(2, 7) )
         HAUTIN    = MOTREA( ADRESS(2, 8) )
         PERDEB    = MOTREA( ADRESS(2, 9) )
         PERFIN    = MOTREA( ADRESS(2,10) )
         PERPAS    = MOTREA( ADRESS(2,11) )
         PERPIC    = MOTREA( ADRESS(2,12) )
         GAMMA     = MOTREA( ADRESS(2,13) )
         TETMIN    = MOTREA( ADRESS(2,14) )
         TETMAX    = MOTREA( ADRESS(2,15) )
         EXPOS     = MOTREA( ADRESS(2,16) )
         RELAX     = MOTREA( ADRESS(2,17) )
         EPSDIS    = MOTREA( ADRESS(2,18) )
         RELDIS    = MOTREA( ADRESS(2,19) )
         ALFABJ    = MOTREA( ADRESS(2,20) )
         GAMMAS    = MOTREA( ADRESS(2,21) )
         KDALLY    = MOTREA( ADRESS(2,22) )
         GDALLY    = MOTREA( ADRESS(2,23) )
         VISCO     = MOTREA( ADRESS(2,24) )
         DIAM90    = MOTREA( ADRESS(2,25) )
         DIAM50    = MOTREA( ADRESS(2,26) )
         MVSED     = MOTREA( ADRESS(2,27) )
         MVEAU     = MOTREA( ADRESS(2,28) )
         FWCOEF    = MOTREA( ADRESS(2,29) )
         RICOEF    = MOTREA( ADRESS(2,30) )
         FFON      = MOTREA( ADRESS(2,31) )
         PMIN      = MOTREA( ADRESS(2,32) )
         PMAX      = MOTREA( ADRESS(2,33) )
C
C MOTS-CLES LOGIQUES :
C
         LISTIN    = MOTLOG( ADRESS(3, 1) )
         INFOGR    = MOTLOG( ADRESS(3, 2) )
         BALAYE    = MOTLOG( ADRESS(3, 3) )
         ALEMON    = MOTLOG( ADRESS(3, 4) )
         ALEMUL    = MOTLOG( ADRESS(3, 5) )
         DEFERL    = MOTLOG( ADRESS(3, 6) )
         FROTTE    = MOTLOG( ADRESS(3, 7) )
         ENTFW     = MOTLOG( ADRESS(3, 8) )
         ENTREG    = MOTLOG( ADRESS(3, 9) )
         ENTRUG    = MOTLOG( ADRESS(3, 10) )
         LISHOU    = MOTLOG( ADRESS(3, 11) )
         VALID     = MOTLOG( ADRESS(3, 12) )
C        EQUATIONS SPHERIQUES CODEES EN DUR
         SPHERI    = .FALSE.
C
C MOTS-CLES CHAINES DE CARACTERES : CERTAINS SONT UTILISES PAR
C                                   LA PROCEDURE DE LANCEMENT
C
         TITCAS    = MOTCAR( ADRESS(4, 1) )(1:72)
         VARDES    = MOTCAR( ADRESS(4, 2) )(1:72)
         CALL MAJUS(VARDES)
         VARIMP    = MOTCAR( ADRESS(4, 3) )(1:72)
         CALL MAJUS(VARIMP)
C        DE 4 A 5 : LU ET UTILISE PAR PRECOS
         ART_FILES(ARTGEO)%NAME=MOTCAR( ADRESS(4,6) )
!        NOMFOR    = MOTCAR( ADRESS(4, 7) )
!        NOMCAS    = MOTCAR( ADRESS(4, 8) )
!         ART_FILES(ARTCAS)%NAME=MOTCAR( ADRESS(4,8) )
         ART_FILES(ARTCLI)%NAME=MOTCAR( ADRESS(4,9) )
!         write(*,*) 'in lecdon ',ART_FILES(ARTGEO)%NAME
         ART_FILES(ARTRES)%NAME=MOTCAR( ADRESS(4,10) )
C        DE 11 A 14 : LU ET UTILISE PAR PRECOS
         ART_FILES(ARTFON)%NAME=MOTCAR( ADRESS(4,15) )
         ART_FILES(ARTBI1)%NAME=MOTCAR( ADRESS(4,16) )
         ART_FILES(ARTBI2)%NAME=MOTCAR( ADRESS(4,17) )
         ART_FILES(ARTFO1)%NAME=MOTCAR( ADRESS(4,18) )
         ART_FILES(ARTFO2)%NAME=MOTCAR( ADRESS(4,19) )
         ART_FILES(ARTRBI)%NAME=MOTCAR( ADRESS(4,20) )
         ART_FILES(ARTRFO)%NAME=MOTCAR( ADRESS(4,21) )
C        DE 22 A 23 : LU ET UTILISE PAR PRECOS OU NON UTILISE
         CDTINI    = MOTCAR( ADRESS(4,24) )(1:72)
         CALL MAJUS(CDTINI)
         BINGEO    = MOTCAR( ADRESS(4,25) )(1:3)
         CALL MAJUS(BINGEO)
         BINRES    = MOTCAR( ADRESS(4,26) )(1:3)
         CALL MAJUS(BINRES)
         ART_FILES(ARTREF)%NAME=MOTCAR( ADRESS(4,28) )
!        FORMAT DU FICHIER DE RESULTATS
         ART_FILES(ARTRES)%FMT = MOTCAR( ADRESS(4,30) )(1:8)
         CALL MAJUS(ART_FILES(ARTRES)%FMT)
!        FORMAT DU FICHIER DE REFERENCE
         ART_FILES(ARTREF)%FMT = MOTCAR( ADRESS(4,31) )(1:8)
         CALL MAJUS(ART_FILES(ARTREF)%FMT)
!        FORMAT DU FICHIER BINAIRE 1
         ART_FILES(ARTBI1)%FMT = MOTCAR( ADRESS(4,32) )(1:8)
         CALL MAJUS(ART_FILES(ARTBI1)%FMT)
!        FORMAT DU FICHIER BINAIRE 2
         ART_FILES(ARTBI2)%FMT = MOTCAR( ADRESS(4,33) )(1:8)
         CALL MAJUS(ART_FILES(ARTBI2)%FMT)
C
      IF(LISTIN) THEN
         IF(LNG.EQ.1) WRITE(LU,101)
         IF(LNG.EQ.2) WRITE(LU,102)
      ENDIF
101   FORMAT(1X,/,19X, '********************************************',/,
     *            19X, '*               LECDON:                    *',/,
     *            19X, '*        APRES APPEL DE DAMOCLES           *',/,
     *            19X, '*     VERIFICATION DES DONNEES LUES        *',/,
     *            19X, '*     SUR LE FICHIER DES PARAMETRES        *',/,
     *            19X, '********************************************',/)
102   FORMAT(1X,/,19X, '********************************************',/,
     *            19X, '*               LECDON:                    *',/,
     *            19X, '*        AFTER CALLING DAMOCLES            *',/,
     *            19X, '*        CHECKING OF DATA  READ            *',/,
     *            19X, '*         IN THE STEERING FILE             *',/,
     *            19X, '********************************************',/)
C
C-----------------------------------------------------------------------
C  NOMS DES VARIABLES POUR LES FICHIERS DE RESULTATS ET DE GEOMETRIE:
C-----------------------------------------------------------------------
C        
C TABLEAU DE LOGIQUES POUR LES SORTIES DE VARIABLES
C        
      CALL NOMVAR_ARTEMIS(TEXTE,TEXTPR,MNEMO)
      CALL SORTIE(VARDES , MNEMO , 100 , SORLEO )
      CALL SORTIE(VARIMP , MNEMO , 100 , SORIMP )
C
C-----------------------------------------------------------------------
C
C DANS LE CAS D'UN BALAYAGE EN PERIODE, LA PREMIERE PERIODE A CALCULER
C EST PERDEB
C
      IF (BALAYE) PER = PERDEB
C
C-----------------------------------------------------------------------
C
C SI IL N'Y A PAS D'IMPRESSIONS SUR LE LISTING, IL NE DOIT PAS Y AVOIR
C D'INFORMATIONS ECRITES SUR LE SOLVEUR
C
      IF (.NOT.LISTIN) INFOGR = .FALSE.
C
C-----------------------------------------------------------------------
C
C POUR UN CALCUL EN HOULE ALEATOIRE, ON INHIBE LES SORTIES GRAPHIQUES ET
C IMPRIMEES DES PHASES, VITESSES, COTE DE LA SURFACE LIBRE, ET
C POTENTIELS CAR ILS NE SIGNIFIENT RIEN.
C
C MAIS ON SORT QUAND MEME UN NOMBRE D'ONDE MOYEN, UNE CELERITE DE
C PHASE MOYENNE, ET UNE VITESSE DE GROUPE MOYENNE CALCULES A PARTIR
C DE LA PERIODE MOYENNE T01, DE MEME QU'UNE INCIDENCE MOYENNE
C
C EN OUTRE, ON VERIFIE QUE LE NOMBRE DE PERIODES ET DE DIRECTIONS
C DE DISCRETISATION EST AU MOINS EGALE A 1
C
      IF(ALEMON .OR. ALEMUL) THEN
C
C       INDICE 2 : PHASE
C
        SORLEO( 2) = .FALSE.
        SORIMP( 2) = .FALSE.
C
C       INDICES 3 ET 4 : U0 ET V0
C
        SORLEO( 3) = .FALSE.
        SORIMP( 3) = .FALSE.
        SORLEO( 4) = .FALSE.
        SORIMP( 4) = .FALSE.
C
C       INDICE 5 : SURFACE LIBRE

C
        SORLEO( 5) = .FALSE.
        SORIMP( 5) = .FALSE.
C
C       INDICES 11 ET 12 : POTENTIEL REEL ET IMAGINAIRE
C
        SORLEO(11) = .FALSE.
        SORIMP(11) = .FALSE.
        SORLEO(12) = .FALSE.
        SORIMP(12) = .FALSE.
C
        IF (NPALE.LE.0) NPALE = 1
        IF (NDALE.LE.0) NDALE = 1
C
      ELSE
C
C       INDICES 17, 18 ET 19 : T01, T02 ET TM
C
        SORLEO(17) = .FALSE.
        SORIMP(17) = .FALSE.
        SORLEO(18) = .FALSE.
        SORIMP(18) = .FALSE.
        SORLEO(19) = .FALSE.
        SORIMP(19) = .FALSE.
      ENDIF
C
C PAS DE BALAYAGE EN PERIODE SI HOULE ALEATOIRE
C
      IF (ALEMON .OR. ALEMUL) BALAYE = .FALSE.
C
C SI HOULE MONODIRECTIONNELLE, IL FAUT NDALE=1
C
      IF (ALEMON .AND. .NOT.ALEMUL) NDALE = 1
C
C-----------------------------------------------------------------------
C
C  PAS DE SORTIE DE VARIABLES PRIVEES SI ELLES N'ONT PAS ETE ALLOUEES
C
      LSTOP = .FALSE.
      DO 15 I=1,4
        IF ((SORLEO(12+I).OR.SORIMP(12+I)).AND.(NPRIV.LT.I)) THEN
          IF (LNG.EQ.1) WRITE(LU,16) I,NPRIV
          IF (LNG.EQ.2) WRITE(LU,17) I,NPRIV
          LSTOP = .TRUE.
        ENDIF
 15   CONTINUE
      IF (LSTOP) STOP
 16   FORMAT(1X,'LA VARIABLE PRIVEE ',1I1,' NE PEUT ETRE UTILISEE '
     *      ,1X,'CAR ELLE N''EST PAS ALLOUEE.',/
     *      ,1X,'AUGMENTER ''NPRIV'' (ACTUELLEMENT ',1I1,' ).',/)
 17   FORMAT(1X,'PRIVATE ARRAY ',1I1,' CANNOT BE USED '
     *      ,1X,'BECAUSE IT WAS NOT ALLOCATED.',/
     *      ,1X,'CHECK ''NPRIV''  (AT THE TIME BEING ',1I1,' ).',/)
C
C-----------------------------------------------------------------------
C
C  IMPRESSION DU TITRE LU
C
      IF(LISTIN) THEN
         IF(LNG.EQ.1) WRITE(LU,3000) TITCAS
         IF(LNG.EQ.2) WRITE(LU,3001) TITCAS
3000     FORMAT(/1X,'TITRE DE L''ETUDE :',1X,A72,/)
3001     FORMAT(/1X,'NAME OF THE STUDY :',1X,A72,/)
      ENDIF
C
C-----------------------------------------------------------------------
C
C PAS DE BANCS DECOUVRANTS ==> MSK = .FALSE.
C
      MSK = .FALSE.
C
C-----------------------------------------------------------------------
C
C OPTION OBLIGATOIRE POUR LE SOLVEUR DIRECT
C
      IF( (SLVART%SLV.EQ.8).AND.(OPTASS.NE.3) ) THEN
          IF(LNG.EQ.1) WRITE(LU,3002)
          IF(LNG.EQ.2) WRITE(LU,3003)
C          
3002    FORMAT(1X,'AVEC SOLVEUR DIRECT, STOCKAGE PAR SEGMENT',/,1X,
     *             'OBLIGATOIRE',///)
3003    FORMAT(1X,'WITH DIRECT SYSTEM SOLVER, EDGE-BASED STORAGE',/,1X,
     *             'IS MANDATORY',///)
          CALL PLANTE(1)
          STOP
      ENDIF    
C
C SOLVEUR DIRECT AVEC PARALLELISME IMPOSSIBLE
C
      IF(NCSIZE.GT.1) THEN
        IF(SLVART%SLV.EQ.8) THEN
           IF(LNG.EQ.1) WRITE(LU,3004)
           IF(LNG.EQ.2) WRITE(LU,3005)
3004       FORMAT(1X,'AVEC PARALLELISME,',/,1X,
     *             'PAS DE SOLVEUR DIRECT',///)
3005       FORMAT(1X,'WITH PARALLELISM,',/,1X,
     *             'NO DIRECT SYSTEM SOLVER',///)
           CALL PLANTE(1)
           STOP
        ENDIF
      ENDIF   
C
C-----------------------------------------------------------------------
C
      RETURN
      END
