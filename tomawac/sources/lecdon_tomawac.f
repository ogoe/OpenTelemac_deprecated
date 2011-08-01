C                       *************************
                        SUBROUTINE LECDON_TOMAWAC
C                       *************************
C
     * (FILE_DESC,PATH,NCAR,CODE)
C
C**********************************************************************
C  TOMAWAC - 6.0    Michel BENOIT (EDF R&D LNHE)            06/12/2004
C**********************************************************************
C
C      FONCTION:
C      =========
C
C    LECTURE DU FICHIER CAS PAR APPEL DU LOGICIEL DAMOCLES.
C
C-----------------------------------------------------------------------
C
C SOUS-PROGRAMME APPELE PAR : HOMERE_TOMAWAC
C SOUS-PROGRAMMES APPELES   : DAMOC 
C
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_TOMAWAC
C
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C DEGRAD = parametre pour la conversion des degres en radians
      DOUBLE PRECISION DEGRAD, PIS2
      PARAMETER(DEGRAD=0.01745329252D0, PIS2=1.570796327D0)
      CHARACTER*8      MNEMO(MAXVAR)
      INTEGER          K
C
C-----------------------------------------------------------------------
C
C DECLARATION DES TABLEAUX POUR L'APPEL DE DAMOCLES
C
      INTEGER, PARAMETER :: NMAX = 300
C
      INTEGER          ADRESS(4,NMAX),DIMEN(4,NMAX)
      DOUBLE PRECISION MOTREA(NMAX)
      INTEGER          MOTINT(NMAX)
      LOGICAL          MOTLOG(NMAX)
      CHARACTER*144    MOTCAR(NMAX)
      CHARACTER*72     MOTCLE(4,NMAX,2)
      INTEGER          TROUVE(4,NMAX)
      LOGICAL          DOC
      CHARACTER(LEN=250) :: NOM_CAS
      CHARACTER(LEN=250) :: NOM_DIC
!ARGUMENTS
      CHARACTER(LEN=24), INTENT(IN)     :: CODE
      CHARACTER(LEN=144), INTENT(INOUT) :: FILE_DESC(4,NMAX)
      INTEGER, INTENT(IN)               :: NCAR
      CHARACTER(LEN=250), INTENT(IN)    :: PATH
      INTEGER :: I
C
C FIN DES VARIABLES AJOUTEES POUR LE NOUVEAU DAMOCLES
C
C
C***********************************************************************
C
      IF (LNG.EQ.1) WRITE(LU,1)
      IF (LNG.EQ.2) WRITE(LU,2)
1     FORMAT(1X,/,19X, '********************************************',/,
     *            19X, '*     SOUS-PROGRAMME LECDON_TOMAWAC        *',/,
     *            19X, '*           APPEL DE DAMOCLES              *',/,
     *            19X, '*     VERIFICATION DES DONNEES LUES        *',/,
     *            19X, '*           SUR LE FICHIER CAS             *',/,
     *            19X, '********************************************',/)
2     FORMAT(1X,/,19X, '********************************************',/,
     *            19X, '*        SUBROUTINE LECDON_TOMAWAC         *',/,
     *            19X, '*           CALL OF DAMOCLES               *',/,
     *            19X, '*        VERIFICATION OF READ DATA         *',/,
     *            19X, '*            ON STEERING FILE              *',/,
     *            19X, '********************************************',/)
C
C-----------------------------------------------------------------------
C
C INITIALISATIONS :
C
      DO K=1,NMAX
C       UN FICHIER NON DONNE PAR DAMOCLES SERA RECONNU PAR UN BLANC
C       (IL N'EST PAS SUR QUE TOUS LES COMPILATEURS INITIALISENT AINSI)
        MOTCAR(K)(1:1)=' '
C
        DIMEN(1,K) = 0
        DIMEN(2,K) = 0
        DIMEN(3,K) = 0
        DIMEN(4,K) = 0
      ENDDO
C
C     IMPRESSION DE LA DOC
      DOC = .FALSE.
C
C-----------------------------------------------------------------------
C     OUVERTURE DES FICHIERS DICTIONNAIRE ET CAS
C-----------------------------------------------------------------------
C
      IF(NCAR.GT.0) THEN
C
        NOM_DIC=PATH(1:NCAR)//'WACDICO'
        NOM_CAS=PATH(1:NCAR)//'WACCAS'
C
      ELSE
C
        NOM_DIC='WACDICO'
        NOM_CAS='WACCAS'
C
      ENDIF
C
      OPEN(2,FILE=NOM_DIC,FORM='FORMATTED',ACTION='READ')
      OPEN(3,FILE=NOM_CAS,FORM='FORMATTED',ACTION='READ')
C
      CALL DAMOCLE
     *( ADRESS, DIMEN , NMAX  , DOC    , LNG   , LU    , MOTINT,
     *  MOTREA, MOTLOG, MOTCAR, MOTCLE , TROUVE, 2  , 3  ,
     *  .FALSE.,FILE_DESC)
C
C     DECRYPTAGE DES CHAINES SUBMIT
C
      CALL READ_SUBMIT(WAC_FILES,44,CODE,FILE_DESC,300)
C
C-----------------------------------------------------------------------
C
C     RETRIEVING FILES NUMBERS IN TOMAWAC FORTRAN PARAMETERS
C     AT THIS LEVEL LOGICAL UNITS ARE EQUAL TO THE FILE NUMBER
C
      DO I=1,44
        IF(WAC_FILES(I)%TELNAME.EQ.'WACGEO') THEN
          WACGEO=I
        ELSEIF(WAC_FILES(I)%TELNAME.EQ.'WACCAS') THEN
          WACCAS=I
        ELSEIF(WAC_FILES(I)%TELNAME.EQ.'WACCLI') THEN
          WACCLI=I
        ELSEIF(WAC_FILES(I)%TELNAME.EQ.'WACFON') THEN
          WACFON=I
        ELSEIF(WAC_FILES(I)%TELNAME.EQ.'WACRES') THEN
          WACRES=I
        ELSEIF(WAC_FILES(I)%TELNAME.EQ.'WACREF') THEN
          WACREF=I
        ELSEIF(WAC_FILES(I)%TELNAME.EQ.'WACLEO') THEN
          WACLEO=I
        ELSEIF(WAC_FILES(I)%TELNAME.EQ.'WACPRE') THEN
          WACPRE=I
        ELSEIF(WAC_FILES(I)%TELNAME.EQ.'WACRBI') THEN
          WACRBI=I
        ELSEIF(WAC_FILES(I)%TELNAME.EQ.'WACCOB') THEN
          WACCOB=I
        ELSEIF(WAC_FILES(I)%TELNAME.EQ.'WACCOF') THEN
          WACCOF=I
        ELSEIF(WAC_FILES(I)%TELNAME.EQ.'WACBI1') THEN
          WACBI1=I
        ELSEIF(WAC_FILES(I)%TELNAME.EQ.'WACFO1') THEN
          WACFO1=I
        ELSEIF(WAC_FILES(I)%TELNAME.EQ.'WACVEB') THEN
          WACVEB=I
        ELSEIF(WAC_FILES(I)%TELNAME.EQ.'WACVEF') THEN
          WACVEF=I
        ELSEIF(WAC_FILES(I)%TELNAME.EQ.'WACPAR') THEN
          WACPAR=I
        ELSEIF(WAC_FILES(I)%TELNAME.EQ.'WACMAB') THEN
          WACMAB=I
        ELSEIF(WAC_FILES(I)%TELNAME.EQ.'WACMAF') THEN
          WACMAF=I
        ENDIF
      ENDDO
C
C-----------------------------------------------------------------------
C
C   AFFECTATION DES PARAMETRES SOUS LEUR NOM EN FORTRAN
C
C-----------------------------------------------------------------------
C
C MOTS CLES ENTIERS
C
      GRAPRD = MOTINT( ADRESS(1,  1) )
      LISPRD = MOTINT( ADRESS(1,  2) )
      NIT    = MOTINT( ADRESS(1,  3) )
      NPLAN  = MOTINT( ADRESS(1,  4) )
      NF     = MOTINT( ADRESS(1,  5) )
      GRADEB = MOTINT( ADRESS(1,  6) )
      LISFON = MOTINT( ADRESS(1,  7) )
      SVENT  = MOTINT( ADRESS(1,  8) )
      SMOUT  = MOTINT( ADRESS(1,  9) )
      SFROT  = MOTINT( ADRESS(1, 10) )
      STRIF  = MOTINT( ADRESS(1, 11) )
      INDIC  = MOTINT( ADRESS(1, 12) )
      INDIV  = MOTINT( ADRESS(1, 13) )
      NSITS  = MOTINT( ADRESS(1, 14) )
      INISPE = MOTINT( ADRESS(1, 15) )
      IDTEL  = MOTINT( ADRESS(1, 16) )
      NPTT   = MOTINT( ADRESS(1, 17) )
      LVMAC  = MOTINT( ADRESS(1, 18) )
      SBREK  = MOTINT( ADRESS(1, 19) )
      IQBBJ  = MOTINT( ADRESS(1, 20) )
      IHMBJ  = MOTINT( ADRESS(1, 21) )
      IFRBJ  = MOTINT( ADRESS(1, 22) )
      IWHTG  = MOTINT( ADRESS(1, 23) )
      IFRTG  = MOTINT( ADRESS(1, 24) )
      IDISRO = MOTINT( ADRESS(1, 25) )
      IEXPRO = MOTINT( ADRESS(1, 26) )
      IFRRO  = MOTINT( ADRESS(1, 27) )
      IFRIH  = MOTINT( ADRESS(1, 28) )
      NDTBRK = MOTINT( ADRESS(1, 29) )
      LIMIT  = MOTINT( ADRESS(1, 30) )
C
C         NPRIV     = MOTINT( ADRESS(1, 32?) ) 'ajouter dans le dico'
C
      STRIA  = MOTINT( ADRESS(1, 32) )
      LIMSPE = MOTINT( ADRESS(1, 33) )
      LAM    = MOTINT( ADRESS(1, 34) )
      INDIM  = MOTINT( ADRESS(1, 35) )
      IDHMA  = MOTINT( ADRESS(1, 36) )
      FRABI  = MOTINT( ADRESS(1, 37) )
      NPRIV  = MOTINT( ADRESS(1, 38) )
      FRABL  = MOTINT( ADRESS(1, 39) )
C     COORDONNEES DU POINT ORIGINE EN (X,Y)
      I_ORIG = MOTINT( ADRESS(1, 40) )
      J_ORIG = MOTINT( ADRESS(1, 40)+1 )
C     MOT-CLEF DEBUGGER
      DEBUG  = MOTINT( ADRESS(1, 41) )
C
C     STANDARD DU FICHIER GEOMETRIE
      STDGEO = 3
C
C MOTS CLES REELS
C
      DT     = MOTREA( ADRESS(2,  1) )
      F1     = MOTREA( ADRESS(2,  2) )
      RAISF  = MOTREA( ADRESS(2,  3) )
      NPLEO  = DIMEN(2,4)
      DO K=1,DIMEN(2,4)
        XLEO(K)= MOTREA( ADRESS(2,  4) + K-1)
      ENDDO
      DO K=1,DIMEN(2,5)
        YLEO(K)= MOTREA( ADRESS(2,  5) + K-1)
      ENDDO
      DDC    = MOTREA( ADRESS(2,  6) )
      CFROT1 = MOTREA( ADRESS(2,  7) )
      CMOUT1 = MOTREA( ADRESS(2,  8) )
      CMOUT2 = MOTREA( ADRESS(2,  9) )
      ROAIR  = MOTREA( ADRESS(2, 10) )
      ROEAU  = MOTREA( ADRESS(2, 11) )
      BETAM  = MOTREA( ADRESS(2, 12) )
      XKAPPA = MOTREA( ADRESS(2, 13) )
      ALPHA  = MOTREA( ADRESS(2, 14) )
      DECAL  = MOTREA( ADRESS(2, 15) )
      ZVENT  = MOTREA( ADRESS(2, 16) )
      CDRAG  = MOTREA( ADRESS(2, 17) )
      HM0    = MOTREA( ADRESS(2, 18) )
      FPIC   = MOTREA( ADRESS(2, 19) )
      GAMMA  = MOTREA( ADRESS(2, 20) )
      SIGMAA = MOTREA( ADRESS(2, 21) )
      SIGMAB = MOTREA( ADRESS(2, 22) )
      ALPHIL = MOTREA( ADRESS(2, 23) )
      FETCH  = MOTREA( ADRESS(2, 24) )
      FREMAX = MOTREA( ADRESS(2, 25) )
      TETA1  = MOTREA( ADRESS(2, 26) )*DEGRAD
      SPRED1 = MOTREA( ADRESS(2, 27) )
      TETA2  = MOTREA( ADRESS(2, 28) )*DEGRAD
      SPRED2 = MOTREA( ADRESS(2, 29) )
      XLAMDA = MOTREA( ADRESS(2, 30) )
      TAILF  = MOTREA( ADRESS(2, 31) )
      E2FMIN = MOTREA( ADRESS(2, 32) )
      ALFABJ = MOTREA( ADRESS(2, 33) )
      GAMBJ1 = MOTREA( ADRESS(2, 34) )
      GAMBJ2 = MOTREA( ADRESS(2, 35) )
      BORETG = MOTREA( ADRESS(2, 36) )
      GAMATG = MOTREA( ADRESS(2, 37) )
      ALFARO = MOTREA( ADRESS(2, 38) )
      GAMARO = MOTREA( ADRESS(2, 39) )
      GAM2RO = MOTREA( ADRESS(2, 40) )
      BETAIH = MOTREA( ADRESS(2, 41) )
      EM2SIH = MOTREA( ADRESS(2, 42) )
      COEFHS = MOTREA( ADRESS(2, 43) )
      XDTBRK = MOTREA( ADRESS(2, 44) )
      XLAMD  = MOTREA( ADRESS(2, 45) )
      ZREPOS = MOTREA( ADRESS(2, 46) )
      ALFLTA = MOTREA( ADRESS(2, 47) )
      RFMLTA = MOTREA( ADRESS(2, 48) )
      KSPB   = MOTREA( ADRESS(2, 49) )
      BDISPB = MOTREA( ADRESS(2, 50) )*DEGRAD
      BDSSPB = MOTREA( ADRESS(2, 51) )*DEGRAD
      HM0L   = MOTREA( ADRESS(2, 52) )
      FPICL  = MOTREA( ADRESS(2, 53) )
      SIGMAL = MOTREA( ADRESS(2, 54) )
      SIGMBL = MOTREA( ADRESS(2, 55) )
      APHILL = MOTREA( ADRESS(2, 56) )
      FETCHL = MOTREA( ADRESS(2, 57) )
      FPMAXL = MOTREA( ADRESS(2, 58) )
      TETA1L = MOTREA( ADRESS(2, 59) )*DEGRAD
      SPRE1L = MOTREA( ADRESS(2, 60) )
      TETA2L = MOTREA( ADRESS(2, 61) )*DEGRAD
      SPRE2L = MOTREA( ADRESS(2, 62) )
      XLAMDL = MOTREA( ADRESS(2, 63) )
      GAMMAL = MOTREA( ADRESS(2, 64) )
      PROMIN = MOTREA( ADRESS(2, 65) )
      VX_CTE = MOTREA( ADRESS(2, 66) )
      VY_CTE = MOTREA( ADRESS(2, 67) )
      CIMPLI = MOTREA( ADRESS(2, 68) )
C
C MOTS CLES LOGIQUES
C
      TSOU   = MOTLOG( ADRESS(3,  1) )
      SPHE   = MOTLOG( ADRESS(3,  2) )
      GLOB   = MOTLOG( ADRESS(3,  3) )
      SUIT   = MOTLOG( ADRESS(3,  4) )
      PROINF = MOTLOG( ADRESS(3,  5) )
      COUSTA = MOTLOG( ADRESS(3,  6) )
      VENT   = MOTLOG( ADRESS(3,  7) )
      DONTEL = MOTLOG( ADRESS(3,  8) )
      PROP   = MOTLOG( ADRESS(3,  9) )
      VENSTA = MOTLOG( ADRESS(3, 10) )
      VALID  = MOTLOG( ADRESS(3, 11) )
      MAREE  = MOTLOG( ADRESS(3, 12) )
      TRIGO  = MOTLOG( ADRESS(3, 13) )
      SPEULI = MOTLOG( ADRESS(3, 14) )
C
C MOTS CLES CHAINES DE CARACTERES
C
      TITCAS = MOTCAR( ADRESS(4, 1) ) (1:72)
      SORT2D = MOTCAR( ADRESS(4, 2) ) (1:72)
C
C FICHIERS ECRITS DANS LE FICHIER CAS
C
      WAC_FILES(WACGEO)%NAME=MOTCAR( ADRESS(4,3) )
!     NOMFOR = MOTCAR( ADRESS(4, 4) )
!     NOMCAS = MOTCAR( ADRESS(4, 5) )
      WAC_FILES(WACCLI)%NAME=MOTCAR( ADRESS(4,6) )
      WAC_FILES(WACFON)%NAME=MOTCAR( ADRESS(4,7) )
      WAC_FILES(WACRES)%NAME=MOTCAR( ADRESS(4,8) )
      WAC_FILES(WACLEO)%NAME=MOTCAR( ADRESS(4,9) )
      WAC_FILES(WACPRE)%NAME=MOTCAR( ADRESS(4,10) )
      WAC_FILES(WACRBI)%NAME=MOTCAR( ADRESS(4,11) )
      WAC_FILES(WACCOB)%NAME=MOTCAR( ADRESS(4,12) )
      WAC_FILES(WACCOF)%NAME=MOTCAR( ADRESS(4,13) )
      WAC_FILES(WACBI1)%NAME=MOTCAR( ADRESS(4,14) )
      WAC_FILES(WACFO1)%NAME=MOTCAR( ADRESS(4,15) )
      BINGEO = MOTCAR( ADRESS(4,16) )(1:3)
      BINRES = MOTCAR( ADRESS(4,17) )(1:3)
      BINLEO = MOTCAR( ADRESS(4,18) )(1:3)
      BINCOU = MOTCAR( ADRESS(4,19) )(1:3)
      BINRBI = MOTCAR( ADRESS(4,20) )(1:3)
      BINPRE = MOTCAR( ADRESS(4,21) )(1:3)
      VERS   = MOTCAR( ADRESS(4,22) )(1:4)
C
C     DE 23 a 28 MOTS CLES POUR LE CRAY NON UTILES ICI
C
      BINVEN = MOTCAR( ADRESS(4,29) )(1:3)
      BINBI1 = MOTCAR( ADRESS(4,30) )(1:3)
      WAC_FILES(WACVEB)%NAME=MOTCAR( ADRESS(4,31) )
      WAC_FILES(WACVEF)%NAME=MOTCAR( ADRESS(4,32) )
      WAC_FILES(WACPAR)%NAME=MOTCAR( ADRESS(4,33) )
      WAC_FILES(WACREF)%NAME=MOTCAR( ADRESS(4,34) )
      WAC_FILES(WACMAB)%NAME=MOTCAR( ADRESS(4,35) )
      WAC_FILES(WACMAF)%NAME=MOTCAR( ADRESS(4,36) )
      BINMAR = MOTCAR( ADRESS(4,37) )(1:3)
      EQUA   = 'TOMAWAC-COWADIS'
!BD_INCKA format d'ecriture des fichiers
!     FORMAT DU FICHIER DE RESULTATS
      WAC_FILES(WACRES)%FMT = MOTCAR( ADRESS(4,40) )(1:8)
      CALL MAJUS(WAC_FILES(WACRES)%FMT)
!     FORMAT DU FICHIER DE RESULTATS DU CALCUL PRECEDENT
!     SEDIMENT...
      WAC_FILES(WACPRE)%FMT = MOTCAR( ADRESS(4,41) )(1:8)
      CALL MAJUS(WAC_FILES(WACPRE)%FMT)
!     FORMAT DU FICHIER DE REFERENCE
      WAC_FILES(WACREF)%FMT = MOTCAR( ADRESS(4,42) )(1:8)
      CALL MAJUS(WAC_FILES(WACREF)%FMT)
!     FORMAT DU FICHIER BINAIRE 1
      WAC_FILES(WACBI1)%FMT = MOTCAR( ADRESS(4,43) )(1:8)
      CALL MAJUS(WAC_FILES(WACBI1)%FMT)
!     FORMAT DU FICHIER DE SPECTRE
      WAC_FILES(WACLEO)%FMT = MOTCAR( ADRESS(4,44) )(1:8)
      CALL MAJUS(WAC_FILES(WACLEO)%FMT)
C
C  CORRECTIONS ET CALCULS D'AUTRES PARAMETRES QUI SE DEDUISENT
C  DE CEUX QUE L'ON VIENT DE LIRE
C
      IF(COUSTA.OR.MAREE) THEN
        COURAN=.TRUE.
      ELSE
        COURAN=.FALSE.
      ENDIF
      IF(.NOT.VENT.AND.SVENT.NE.0) THEN
         IF(LNG.EQ.1) 
     *         WRITE(LU,*)
     *      'INCOHERENCE DES MOTS CLES DU VENT => PAS DE VENT'
         IF(LNG.EQ.2) 
     *         WRITE(LU,*)
     *      'INCOMPATIBILITY OF KEY-WORDS CONCERNING WIND => NO
     * WIND'
        SVENT=0
      ENDIF
      IF(TRIGO) THEN
        TETA1  = PIS2-TETA1
        TETA2  = PIS2-TETA2
        TETA1L = PIS2-TETA1L
        TETA2L = PIS2-TETA2L
        BDISPB = PIS2-BDISPB
        BDSSPB = PIS2-BDSSPB
      ENDIF
      IF(LIMIT.LT.0.OR.LIMIT.GT.2) THEN
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'INCOHERENCE DE L OPTION DE LIMITEUR'
          WRITE(LU,*) 'VALEUR LUE = ',LIMIT
          WRITE(LU,*) 'ON PREND LA VALEUR PAR DEFAUT LIMIT=1'
        ENDIF  
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) 'INCOMPATIBILITY OF LIMITER OPTION'
          WRITE(LU,*) 'VALUE READ = ',LIMIT
          WRITE(LU,*) 'WE TAKE THE DEFAULT VALUE LIMIT=1'
        ENDIF
        LIMIT=1
      ENDIF
      IF(CIMPLI.LT.0.OR.CIMPLI.GT.1) THEN
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'INCOHERENCE DU COEFFICIENT D IMPLICITATION'
          WRITE(LU,*) 'VALEUR LUE = ',CIMPLI
          WRITE(LU,*) 'ON PREND LA VALEUR PAR DEFAUT CIMPLI=0.5'
        ENDIF  
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) 'INCOMPATIBILITY OF IMPLICITATION COEFFICIENT'
          WRITE(LU,*) 'VALUE READ = ',CIMPLI
          WRITE(LU,*) 'WE TAKE THE DEFAULT VALUE CIMPLI=0.5'
        ENDIF
        CIMPLI=0.5D0
      ENDIF     
C
C
C-----------------------------------------------------------------------
C  NOMS DES VARIABLES POUR LES FICHIERS DE RESULTATS ET DE GEOMETRIE:
C-----------------------------------------------------------------------
C
C TABLEAU DE LOGIQUES POUR LES SORTIES DE VARIABLES
C
      CALL NOMVAR_TOMAWAC(TEXTE,TEXTPR,MNEMO,MAXVAR)
C      
C$DC$ BUG : tableaux MNEMO et SORLEO DIMENionnes a MAXVAR      
C             tres inferieur a 100 !          
      CALL SORTIE(SORT2D , MNEMO , MAXVAR , SORLEO )
C
C.....Si pas de vent, on n'imprime pas les composantes du vent en sortie.
      IF (.NOT.VENT) THEN
        SORLEO( 9)=.FALSE.
        SORLEO(10)=.FALSE.
      ENDIF
C
C.....Si profondeur infinie, pas de calcul des contraintes de radiation.
      IF (PROINF) THEN
        IF (SORLEO(11) .OR. SORLEO(12) .OR. SORLEO(13) .OR.
     *      SORLEO(14) .OR. SORLEO(15) ) THEN
           IF (LNG.EQ.1) THEN
             WRITE(LU,*) '******************************************'
             WRITE(LU,*) ' LE CALCUL DES CONTRAINTES DE RADIATION ET'
             WRITE(LU,*) '  DES FORCES MOTRICES NE S''EFFECTUE PAS'
             WRITE(LU,*) '      PAS EN PROFONDEUR INFINIE          '
             WRITE(LU,*) '******************************************'
           ELSE
             WRITE(LU,*) '*****************************************'
             WRITE(LU,*) '   RADIATION STRESSES ARE NOT COMPUTED  '
             WRITE(LU,*) '       OVER INFINITE WATER DEPTHS       '
             WRITE(LU,*) '******************************************'
           ENDIF
           DO K=11,15
             SORLEO(K) = .FALSE.
           ENDDO         
        ENDIF
      ENDIF
C
      DO K=1,MAXVAR
        SORIMP(K)=.FALSE.
      ENDDO
C
C
 1001 FORMAT('*** INCOMPATIBILITE DES MOTS CLES ***
     *                ARRET DU PROGRAMME')
 1002 FORMAT('*** INCOMPATIBILITY OF THE KEY WORDS ***
     *                   PROGRAM STOP')
C
C-----------------------------------------------------------------------
C
      RETURN
      END
