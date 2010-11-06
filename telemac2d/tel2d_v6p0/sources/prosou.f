C                       *****************
                        SUBROUTINE PROSOU
C                       *****************
C
     *(FU,FV,SMH,    UN,VN,HN,GRAV,NORD,
     * FAIR,WINDX,WINDY,VENT,HWIND,CORIOL,FCOR,
     * SPHERI,YASMH,COSLAT,SINLAT,AT,LT,
     * NREJET,NREJEU,DSCE,ISCE,T1,MESH,MSK,MASKEL,
     * MAREE,MARDAT,MARTIM,PHI0,OPTSOU,COUROU,NPTH,VARCL,NVARCL,VARCLA,
     * UNSV2D)
C
C***********************************************************************
C  TELEMAC 2D VERSION 6.0   08/04/08  J-M HERVOUET (LNHE) 01 30 87 80 18
C                                          
C***********************************************************************
C
C  FONCTION :
C
C     PREPARATION DES TERMES SOURCES :
C
C     DANS L'EQUATION DE CONTINUITE
C
C     DANS LES EQUATIONS DYNAMIQUES
C
C
C     SONT PRIS EN COMPTE : - LE VENT.
C                           - LA FORCE DE CORIOLIS.
C                           - LA FORCE GENERATRICE DE LA MAREE.
C                           - LES SOURCES ET LES PUITS.
C     NOTA:
C
C     LE FROTTEMENT DE FOND EST PRIS EN COMPTE DANS LA PROPAGATION
C     PAR L'APPEL A FROTXY, IL EST SEMI-IMPLICITE.
C
C     SI L'ON AJOUTE DES TERMES SOURCES OU PUITS A L'EQUATION DE
C     CONTINUITE, IL FAUT LE DIRE EN METTANT LA VARIABLE YASMH A
C     TRUE
C
C     ON CONSTRUIT D'ABORD DES TERMES SOURCES FU ET FV EN P1
C     PUIS ON LES ETEND AU QUASI-BULLE SI IL LE FAUT.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   FU           |<-- | TERMES SOURCE TRAITE EN P1 SUR L'EQUATION EN U
C |   FV           |<-- | TERMES SOURCE TRAITE EN P1 SUR L'EQUATION EN V
C |   SMH          |<-- | TERME SOURCE DANS L'EQUATION DE CONTINUITE
C |   UN , VN      | -->| COMPOSANTES DES VECTEURS VITESSES A TN.
C |   HN           | -->| HAUTEURS D'EAU A TN .
C |   GRAV         | -->| PESANTEUR.
C |   NORD         | -->| DIRECTION DU NORD EN DEGRES PAR RAPPORT A
C |                |    | L'AXE DES Y (SENS TRIGONOMETRIQUE)
C |   FAIR         | -->| COEFFICIENT DE FROTTEMENT DE L'AIR.
C |   WINDX,Y      | -->| VITESSE DU VENT EN SURFACE.
C |   VENT         | -->| PRISE EN COMPTE DES EFFORTS DUS AU VENT .
C |   HWIND        | -->| MINIMUM DEPTH FOR TAKING WIND INTO ACCOUNT
C |   CORIOL       | -->| PRISE EN COMPTE DES EFFORTS DE CORIOLIS .
C |   FCOR         | -->| PARAMETRE DE CORIOLIS.
C |   SPHERI       | -->| =TRUE : COORDONNEES SPHERIQUES
C |   YASMH        |<-->| =TRUE SI SMH NON NUL. AJOUTE UN TERME SOURCE
C |                |    | IMPLICITE DANS L'EQUATION DU TRACEUR
C |   COSLAT       | -->| COS DE LA LATITUDE EN COORDONNEES SPHERIQUES
C |   SINLAT       | -->| SIN DE LA LATITUDE EN COORDONNEES SPHERIQUES
C |   AT           | -->| TIME
C |   LT           | -->| TIME STEP NUMBER
C |   NREJET       | -->| NOMBRE DE POINTS SOURCES
C |   NREJEU       | -->| NOMBRE DE VITESSES DES SOURCES DONNEES
C |                |    | SI NREJEU=0 ON CONSIDERE QUE LA VITESSE DES
C |                |    | SOURCES EST EGALE A CELLE DU COURANT.
C |   ISCE,DSCE    | -->| POINTS SOURCES, DEBITS DE LA SOURCE
C |   T1           | -->| TABLEAU DE TRAVAIL
C |   MESH         | -->| MAILLAGE
C |   MSK          | -->|  SI OUI, PRESENCE D'ELEMENTS MASQUES.
C |   MASKEL       | -->|  TABLEAU DE MASQUAGE DES ELEMENTS
C |                |    |  =1. : NORMAL   =0. : ELEMENT MASQUE
C |   COUROU       | -->| IF YES, WAVE DRIVEN CURRENTS TAKEN INTO ACCOUNT
C |   NPTH         | -->| RECORD NUMBER IN THE WAVE CURRENTS FILE
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELMAC
C
C SOUS-PROGRAMME APPELE : OV
C
C***********************************************************************
C
C    LES TERMES RESPECTIFS SONT:
C    ==========================
C
C
C     * LE VENT .
C       ---------
C
C                                 1                         2      2
C                FU           =  --- * F    * U    * SQRT( U    + V   )
C                  VENT           H     AIR    AIR          AIR    AIR
C
C                                 1                         2      2
C                FV           =  --- * F    * V    * SQRT( U    + V   )
C                  VENT           H     AIR    AIR          AIR    AIR
C
C           AVEC :
C
C                  UAIR   :  VITESSE DU VENT DANS LA DIRECTION X
C                  VAIR   :  VITESSE DU VENT DANS LA DIRECTION Y
C                  FAIR   :  COEFFICIENT DE FROTTEMENT DE L'AIR
C
C
C     * LA FORCE DE CORIOLIS.
C       ---------------------
C
C                FU           =  + FCOR * V
C                  CORIOLIS
C
C                FV           =  - FCOR * U
C                  CORIOLIS
C
C           AVEC :
C
C                  U       :  VITESSE DU FLUIDE DANS LA DIRECTION X
C                  V       :  VITESSE DU FLUIDE DANS LA DIRECTION Y
C                  FCOR    :  PARAMETRE DE CORIOLIS
C
C
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_TELEMAC2D, ONLY : T2D_FILES,T2DBI1
      USE INTERFACE_TELEMAC2D, EX_PROSOU => PROSOU
! --- JP Renaud start ---
      USE M_COUPLING_ESTEL3D
! --- JP Renaud end ---
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C     TABLEAUX DE TRAVAIL
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: T1
C
C-----------------------------------------------------------------------
C
C     VECTEURS
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: FU,FV,SMH
      TYPE(BIEF_OBJ), INTENT(IN)    :: MASKEL,UN,VN,HN,UNSV2D
      TYPE(BIEF_OBJ), INTENT(IN)    :: WINDX,WINDY,COSLAT,SINLAT
C
C-----------------------------------------------------------------------
C
C     STRUCTURE DE MAILLAGE
C
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
C
C-----------------------------------------------------------------------
C
      INTEGER, INTENT(IN)           :: NVARCL,LT,NREJET,NREJEU,OPTSOU
      INTEGER, INTENT(IN)           :: NPTH
      INTEGER, INTENT(IN)           :: MARDAT(3),MARTIM(3),ISCE(NREJET)
      DOUBLE PRECISION, INTENT(IN)  :: HWIND,AT,FAIR,FCOR,DSCE(NREJET)
      DOUBLE PRECISION, INTENT(IN)  :: GRAV,NORD,PHI0
      CHARACTER(LEN=32), INTENT(IN) :: VARCLA(NVARCL)
      LOGICAL, INTENT(IN)           :: VENT,MAREE,CORIOL,SPHERI,MSK
      LOGICAL, INTENT(IN)           :: COUROU
      LOGICAL, INTENT(INOUT)        :: YASMH
      TYPE(BIEF_OBJ), INTENT(INOUT) :: VARCL
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER N,I,IELMU,IELMH,IELM1,NPOIN,IR,ERR,NP
C
      DOUBLE PRECISION PI,WROT,WD,ATH
C      
      CHARACTER*16 NOMX,NOMY
      LOGICAL DEJALU,OKX,OKY,OKC
      DATA DEJALU /.FALSE./
      REAL, ALLOCATABLE :: W(:)
      TYPE(BIEF_OBJ) :: FXH,FYH
      SAVE FXH,FYH,W
C
      INTRINSIC SQRT,MAX,ACOS
C
C-----------------------------------------------------------------------
C  EXTRACTION DES COORDONNEES X, DU NOMBRE DE POINTS P1
C                                ET DE L'ELEMENT P1 DU MAILLAGE.
C-----------------------------------------------------------------------
C
      IELM1 = MESH%X%ELM
      NPOIN = MESH%NPOIN
C
C-----------------------------------------------------------------------
C  INITIALISATION
C-----------------------------------------------------------------------
C
      CALL CPSTVC(UN,FU)
      CALL CPSTVC(VN,FV)
      CALL OS( 'X=0     ' , X=FU )
      CALL OS( 'X=0     ' , X=FV )
C
C-----------------------------------------------------------------------
C
C  CALCUL AVEC VENT
C  ----------------
C
C                               1                         2     2
C              FU           =  --- * F    * U    * SQRT( U   + V    )
C                VENT           H     AIR    AIR          AIR   AIR
C
C
C                               1                         2     2
C              FV           =  --- * F    * V    * SQRT( U   + V    )
C                VENT           H     AIR    AIR          AIR   AIR
C
C
      IF(VENT) THEN
C
C  TRAITEMENT PROVISOIRE DES BANCS DECOUVRANTS
C  ON NE PREND EN COMPTE L'EFFET DU VENT QUE SI LA PROFONDEUR EST
C  SUPERIEURE A 1 M.
C
C  ICI ON SUPPOSE QUE LE VENT EST DONNE EN P1.
C
        DO 10 N=1,NPOIN
          IF (HN%R(N).GT.HWIND) THEN
            WD = SQRT( WINDX%R(N)**2 + WINDY%R(N)**2 )
            FU%R(N) = FU%R(N) + FAIR * WINDX%R(N) * WD / HN%R(N)
            FV%R(N) = FV%R(N) + FAIR * WINDY%R(N) * WD / HN%R(N)
          ENDIF
10      CONTINUE
C
      ENDIF
C
C***********************************************************************
C
C     * AVEC LA FORCE DE CORIOLIS.
C       --------------------------
C
C                FU           =  + FCOR * V
C                  CORIOLIS
C
C                FV           =  - FCOR * U
C                  CORIOLIS
C
      IF(CORIOL) THEN
C
      PI = ACOS(-1.D0)
C
        IF(SPHERI) THEN
C
          WROT = 2 * PI / 86164.D0
          DO 20 I=1,NPOIN
C           FORMULE INDEPENDANTE DE LA DIRECTION DU NORD
            FU%R(I) = FU%R(I) + VN%R(I) * 2 * WROT * SINLAT%R(I)
            FV%R(I) = FV%R(I) - UN%R(I) * 2 * WROT * SINLAT%R(I)
20        CONTINUE
C
C         PRISE EN COMPTE DE LA FORCE GENERATRICE DE LA MAREE
C
          IF(MAREE) THEN
            CALL MARAST(MARDAT,MARTIM,PHI0,NPOIN,AT,
     *                  FU%R,FV%R,MESH%X%R,SINLAT%R,COSLAT%R,GRAV)
          ENDIF
C
          IF(LT.EQ.1) THEN
            IF(LNG.EQ.1) WRITE(LU,11)
            IF(LNG.EQ.2) WRITE(LU,12)
          ENDIF
11        FORMAT(1X,'PROSOU : EN COORDONNEES SHERIQUES, LE',/,
     *           1X,'         COEFFICIENT DE CORIOLIS EST',/,
     *           1X,'         CALCULE EN FONCTION DE LA LATITUDE.',/,
     *           1X,'         LE MOT-CLE ''COEFFICIENT DE CORIOLIS''',/,
     *           1X,'         N''EST DONC PAS PRIS EN COMPTE.')
12        FORMAT(1X,'PROSOU : IN SPHERICAL COORDINATES, THE CORIOLIS',/,
     *           1X,'         PARAMETER DEPENDS ON THE LATITUDE.',/,
     *           1X,'         THE KEY WORD ''CORIOLIS COEFFICIENT''',/,
     *           1X,'         IS CONSEQUENTLY IGNORED.')
C
        ELSE
C
          CALL OS( 'X=X+CY  ' , FU , VN , VN ,  FCOR )
          CALL OS( 'X=X+CY  ' , FV , UN , UN , -FCOR )
C
          IF(LT.EQ.1) THEN
            IF(LNG.EQ.1) WRITE(LU,21)
            IF(LNG.EQ.2) WRITE(LU,22)
          ENDIF
21        FORMAT(1X,'PROSOU : EN COORDONNEES CARTESIENNES, LE',/,
     *           1X,'         COEFFICIENT DE CORIOLIS EST LU DANS LE',/,
     *           1X,'         FICHIER DES PARAMETRES ET CORRESPOND',/,
     *           1X,'         AU MOT-CLE ''COEFFICIENT DE CORIOLIS''',/,
     *           1X,'         IL EST ALORS CONSTANT EN ESPACE')
22        FORMAT(1X,'PROSOU : IN CARTESIAN COORDINATES, THE CORIOLIS',/,
     *           1X,'         PARAMETER IS READ IN THE STEERING FILE',/,
     *           1X,'         IT IS THE KEY WORD ''CORIOLIS',/,
     *           1X,'         COEFFICIENT'', IT IS UNIFORM IN SPACE')
C
        ENDIF
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C  LES SECONDS MEMBRES SONT MIS A LA BONNE DISCRETISATION
C
      IELMU=UN%ELM
C
      IF(IELMU.NE.IELM1) THEN
        CALL CHGDIS(FU,IELM1,IELMU,MESH)
        CALL CHGDIS(FV,IELM1,IELMU,MESH)
      ENDIF
C
C-----------------------------------------------------------------------
C
      IELMH=HN%ELM
      CALL CPSTVC(HN,SMH)
      CALL OS( 'X=0     ' , X=SMH )
C
      YASMH = .FALSE.
C
      IF(NREJET.NE.0) THEN
C
C LE LOGIQUE YASMH PASSE A TRUE
C
      YASMH = .TRUE.
C
C  TERMES SOURCES DANS L'EQUATION DE CONTINUITE
C              ET DANS LA QUANTITE DE MOUVEMENT :
C
C  ATTENTION, SMH EST UTILISE AUSSI POUR LE TRACEUR.
C
C     CALCUL DU VOLUME DES BASES
C     HN EST MIS ICI POUR UNE STRUCTURE BIDON
      CALL VECTOR(T1,'=','MASBAS          ',IELMH,
     *            1.D0,HN,HN,HN,HN,HN,HN,MESH,MSK,MASKEL)
C
      IF(NCSIZE.GT.1) CALL PARCOM(T1,2,MESH)
C
      DO I = 1 , NREJET
C
        IR = ISCE(I)
C       LE TEST SERT POUR LE PARALLELISME, QUAND LE POINT SOURCE
C       N'EST PAS DANS LE SOUS-DOMAINE TRAITE
        IF(IR.GT.0) THEN
         IF(OPTSOU.EQ.1) THEN
C          VERSION "NORMAL"
           SMH%R(IR)=SMH%R(IR)+DSCE(I)/T1%R(IR)
         ELSE
C          VERSION "DIRAC"
           SMH%R(IR) = SMH%R(IR)+DSCE(I)
         ENDIF
        ENDIF
C
      ENDDO
C
C-----------------------------------------------------------------------
C TRAITEMENT EXPLICITE DES APPORTS DE QUANTITE DE MOUVEMENT
C AUX SOURCES
C
      IF(NREJEU.GT.0) THEN
C
      DO I = 1 , NREJEU
C
        IR = ISCE(I)
C       LE TEST SERT POUR LE PARALLELISME, QUAND LE POINT SOURCE
C       N'EST PAS DANS LE SOUS-DOMAINE TRAITE
        IF(IR.GT.0) THEN
C       QUANTITE DE MOUVEMENT APPORTEE PAR LA SOURCE
C      -QUANTITE DE MOUVEMENT PRISE PAR LA SOURCE.
        FU%R(IR)=FU%R(IR) + (VUSCE(AT,I)-UN%R(IR))*
     *  DSCE(I)/(T1%R(IR)*MAX(HN%R(IR),0.1D0))
        FV%R(IR)=FV%R(IR) + (VVSCE(AT,I)-VN%R(IR))*
     *  DSCE(I)/(T1%R(IR)*MAX(HN%R(IR),0.1D0))
        ENDIF
C
      ENDDO
C
      ENDIF
C
      ENDIF
C
C***********************************************************************
C
C     * WITH WAVE DRIVEN CURRENTS
C       -------------------------------------
C
C                FU        =  FXH
C                  COUROU
C
C                FV        =  FYH
C                  COUROU
C
C       FXH ET FYH ARE TAKEN IN A RESULTS FILE FROM 
C       (ARTEMIS, COWADIS, OU TOMAWAC)
C
C       ATTENTION : 1. MESHES MUST BE THE SAME
C       ---------    
C
C                   2. STATIONARY FORCING
C
      IF(COUROU) THEN
C
         IF(.NOT.DEJALU) THEN
C
            ALLOCATE(W(NPOIN),STAT=ERR)
            IF(ERR.NE.0) THEN
              IF(LNG.EQ.1) THEN
                WRITE(LU,*) 'ERREUR D''ALLOCATION DE W DANS PROSOU'
              ENDIF
              IF(LNG.EQ.2) THEN
                WRITE(LU,*) 'MEMORY ALLOCATION ERROR OF W IN PROSOU'
              ENDIF
            ENDIF
            CALL ALLVEC(1,FXH   ,'FXH   ',IELMU, 1 , 2 )
            CALL ALLVEC(1,FYH   ,'FYH   ',IELMU, 1 , 2 )
C
C           NBI1 : FICHIER DE DONNEES BINAIRE 1 ; BINARY DATA FILE 1
            NOMX='FORCE FX        '
            NOMY='FORCE FY        '
            CALL FIND_IN_SEL(FXH,NOMX,T2D_FILES(T2DBI1)%LU,
     *                       W,OKX,NPTH,NP,ATH)
            CALL FIND_IN_SEL(FYH,NOMY,T2D_FILES(T2DBI1)%LU,
     *                       W,OKY,NPTH,NP,ATH)
            IF(.NOT.OKX.OR..NOT.OKY) THEN
C             SECOND TRY (OLD VERSIONS OF ARTEMIS OR TOMAWAC)
              NOMX='FORCE_FX        '
              NOMY='FORCE_FY        '
              CALL FIND_IN_SEL(FXH,NOMX,T2D_FILES(T2DBI1)%LU,
     *                         W,OKX,NPTH,NP,ATH)
              CALL FIND_IN_SEL(FYH,NOMY,T2D_FILES(T2DBI1)%LU,
     *                         W,OKY,NPTH,NP,ATH)
            ENDIF
C           VARIABLES CLANDESTINES TRANSMISES DE TOMAWAC A SISYPHE
            IF(NVARCL.GT.0) THEN
              DO I=1,NVARCL
              CALL FIND_IN_SEL(VARCL%ADR(I)%P,VARCLA(I)(1:16),
     *                         T2D_FILES(T2DBI1)%LU,
     *                         W,OKC,NPTH,NP,ATH)
              IF(.NOT.OKC) THEN
                IF(LNG.EQ.1) WRITE(LU,7) VARCLA(I)(1:16)
                IF(LNG.EQ.2) WRITE(LU,8) VARCLA(I)(1:16)
7             FORMAT(1X,'PROSOU : VARIABLE CLANDESTINE :',/,1X,A16,/,1X,
     *                  '         NON TROUVEE',/,1X,
     *                  '         DANS LE FICHIER DE HOULE')
8             FORMAT(1X,'PROSOU : CLANDESTINE VARIABLE:',/,1X,A16,/,1X,
     *                  '         NOT FOUND',/,1X,
     *                  '         IN THE WAVE RESULTS FILE')
              CALL PLANTE(1)
              STOP
              ENDIF
              ENDDO
            ENDIF
C
            IF(.NOT.OKX.OR..NOT.OKY) THEN
              IF(LNG.EQ.1) WRITE(LU,5)
              IF(LNG.EQ.2) WRITE(LU,6)
 5            FORMAT(1X,'PROSOU : FORCE FX OU FY NON TROUVES',/,1X,
     *                  '         DANS LE FICHIER DE HOULE')
 6            FORMAT(1X,'PROSOU: FORCE FX OR FY NOT FOUND',/,1X,
     *                  '         IN THE WAVE RESULTS FILE')
              CALL PLANTE(1)
              STOP
            ENDIF
            IF(NP.NE.NPOIN) THEN
              IF(LNG.EQ.1) WRITE(LU,95)
              IF(LNG.EQ.2) WRITE(LU,96)
 95           FORMAT(1X,'PROSOU : SIMULATION DES COURANTS DE HOULE.',/,
     *               1X,'LES MAILLAGES HOULE ET COURANTS SONT ',/,
     *               1X,'DIFFERENTS : PAS POSSIBLE POUR LE MOMENT.')
 96           FORMAT(1X,'PROSOU: WAVE DRIVEN CURRENTS MODELLING.',/,
     *               1X,'WAVE AND CURRENT MODELS MESHES ARE ',/,
     *               1X,'DIFFERENT : NOT POSSIBLE AT THE MOMENT.')
C
              CALL PLANTE(1)
              STOP
            ENDIF
C           IMPRESSION SUR LE LISTING
            IF(LNG.EQ.1) WRITE(LU,115) ATH 
            IF(LNG.EQ.2) WRITE(LU,116) ATH
115         FORMAT(1X,/,1X,'PROSOU : COURANTS DE HOULE',/,
     *                  1X,'         LECTURE AU TEMPS ',F10.3,/)
116         FORMAT(1X,/,1X,'PROSOU: WAVE DRIVEN CURRENTS MODELLING',/,
     *                  1X,'         READING FILE AT TIME ',F10.3,/)
            IF(IELMU.NE.IELM1) THEN
              CALL CHGDIS(FXH,IELM1,IELMU,MESH)
              CALL CHGDIS(FYH,IELM1,IELMU,MESH)
            ENDIF
            DEJALU = .TRUE.
C
         ENDIF
C
C        ADDING INTO FU AND FV
C
         CALL OS('X=X+Y   ',X=FU,Y=FXH)
         CALL OS('X=X+Y   ',X=FV,Y=FYH)
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C
C  PRISE EN COMPTE DES INFILTRATIONS DANS LE SOUS-SOL
C  COMMUNICATION AVEC ESTEL-3D
C
C     GET SOURCE TERM FROM ESTEL-3D TO ACCOUNT FROM INFILTRATION
C     CALL TO THE INFILTRATION ROUTINE
C
      CALL INFILTRATION_GET(SMH%R,UNSV2D%R,YASMH)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
