c                       ******************
                        SUBROUTINE ARTEMIS
C                       ******************
C
C***********************************************************************
C
C ARTEMIS VERSION 6.0 21/06/10  PARALLEL VERSION   C. DENIS (SINETICS)             
C 
C ARTEMIS VERSION 5.1 21/04/99  D. AELBRECHT  (LNH) 01.30.87.74.12
C
C  LINKED TO BIEF VERS. 5.0  17/08/94  J-M HERVOUET (LNH) 01 30 87 80 18
C
C***********************************************************************
C
C               AAA  RRRR  TTTTT EEEEE M   M IIIII  SSSS
C              A   A R   R   T   E     MM MM   I   S
C              AAAAA RRRR    T   EEEE  M M M   I    SSS
C              A   A R   R   T   E     M   M   I       S
C              A   A R   R   T   EEEEE M   M IIIII SSSS
C
C         ******************************************************
C         *                                                    *
C         *            RESOLUTION DE L'EQUATION DE             *
C         *                 BERKHOFF MODIFIEE                  *
C         *                                                    *
C         ******************************************************
C
C***********************************************************************
C
C SOUS-PROGRAMME APPELE PAR : HOMERE_ARTEMIS
C
C***********************************************************************
C
C-----------------------------------------------------------------------
C                    DECLARATION DES TYPES ET DIMENSIONS
C-----------------------------------------------------------------------
C
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_ARTEMIS
      USE GRACESTOP
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C ENTIERS
C
      INTEGER LT,NPERBA,I,J
      INTEGER NELBRD,NPFMAX,NELBRX
      INTEGER LPER,LDIR
      INTEGER ALIRE(MAXVAR)
C
C VARIABLE POUR SUB/ROUTINE DISMOY
C
      INTEGER LISHHO
C
C SCALAIRES REELS
C
      DOUBLE PRECISION RADDEG,HIST(1)
C
C VARIABLES POUR LES APPELS DE SOUS-PROGRAMMES DE TELEMAC-2D
C
      INTEGER NVARCL,ISTO
      DOUBLE PRECISION LAMBD0
      LOGICAL RESU,FROVAR,PROLIN,TRAC
C
C POUR PASSAGE D'ARGUMENTS VIDES
C
      INTEGER IBID
      DOUBLE PRECISION BID


      INTEGER  P_IMAX,P_IMIN
      DOUBLE PRECISION P_DMIN
      EXTERNAL P_IMAX,P_IMIN,P_DMIN
      
C
      DATA HIST /9999.D0/
C
C-----------------------------------------------------------------------
C
C  VARIABLES A LIRE EN CAS DE SUITE :
C  0 : ON NE LIT PAS    1 : ON LIT  (VOIR SS-PG NOMVAR)
C
      DATA ALIRE /1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     *            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     *            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     *            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
C
C-----------------------------------------------------------------------
C
      RADDEG = 180.D0/3.141592654D0
C
C=======================================================================
C
C : 1          LECTURE , PREPARATION ET CONTROLE DES DONNEES
C
C=======================================================================
C
C  TYPES DE DISCRETISATION :
C
C  TRIANGLES : P1
      IELM  = 11
C  SEGMENTS  : P1 POUR LA FRONTIERE
      IELMB = 1
C
   
      
C  DIMENSIONS MAXI (CAS D'UN MAILLAGE ADAPTATIF)
C  CES PARAMETRES SONT UTILISES POUR LES APPELS DE BIEF
C
C     POINTS
      NPMAX = NPOIN
C     ELEMENTS
      NELMAX = NELEM
C     ELEMENTS DE BORD
      NELBRD = NPTFR
C     ELEMENTS DE BORD (NOMBRE MAXIMUM)
      NPFMAX = NPTFR
C     POINTS DE BORD
      NELBRX = NPTFR
C
      IF(BALAYE) THEN
        NPERBA = INT((PERFIN-PERDEB)/PERPAS) + 1
      ENDIF
C
C=======================================================================
C
      RESU   = .TRUE.
      FROVAR = .FALSE.
      PROLIN = .FALSE.
      SPHERI = .FALSE.
      TRAC   = .FALSE.
      NVARCL = 0
C
C DANS TELEMAC-2D, SI LIHBOR VAUT KINC, LIHBOR EST FORCE A KSORT
C ON CHANGE DONC POUR PREDA2 LA VALEUR DE KINC AFIN D'EVITER CE FORCAGE
C PAR AILLEURS, DANS TELEMAC-2D, SI LIHBOR VAUT KADH (NON DEFINI ICI) 
C ON IMPRIME UN MESSAGE. POUR ETRE SUR DE L'EVITER, ON PASSE AUSSI ISTO 
C A LA PLACE DE KADH.
C
      ISTO = 100
C
C-----------------------------------------------------------------------
C
C LECTURE DES CONDITIONS LIMITES ET INDICES DES POINTS FRONTIERES.
C
      CALL LECLIM_ARTEMIS
     *(LIHBOR%I,LIUBOR%I,MESH%NPTFR,MESH%NBOR%I,STDGEO,
     * ART_FILES(ARTCLI)%LU,
     * MESH%ISEG%I,MESH%XSEG%R,MESH%YSEG%R,MESH%NACHB%I,NUMLIQ%I,
     * MESH%IFAPAR%I)

C
C-----------------------------------------------------------------------
C
C COMPLEMENT DE LA STRUCTURE DE DONNEES POUR BIEF
C
      CALL INBIEF(LIHBOR%I,KLOG,IT1,IT2,IT3,LVMAC,IELM,
     *         LAMBD0,SPHERI,MESH,T1,T2,OPTASS,PRODUC,EQUA)

      IF (NCSIZE .LE. 1) THEN
         NPOIN_TOT=MESH%NPOIN
         ALLOCATE(XT(NPOIN_TOT))
         ALLOCATE(YT(NPOIN_TOT))
      ENDIF
C-----------------------------------------------------------------------
C  RECHERCHE DU FOND ET DU FROTTEMENT DANS LE FICHIER DE GEOMETRIE :
C-----------------------------------------------------------------------
C     
      CALL FONSTR(T1,ZF,T2,FW,ART_FILES(ARTGEO)%LU,ART_FILES(ARTFON)%LU,
     *            ART_FILES(ARTFON)%NAME,MESH,FFON,LISTIN)
 
C-----------------------------------------------------------------------
C
C PREPARATION DU FICHIER DE RESULTATS (FACULTATIF)
C
C     STANDARD SELAFIN
C
        ! CREATION DU JEU DE DONNEES POUR UN FORMAT DE FICHIER
        ! FORMAT_RES.
        ! LE JEU DE DONNEES EST CREE DANS LE FICHIER NRES, ET EST
        ! DEFINIT PAR UN TITRE ET LES VARIABLES A ECRIRE. 
        CALL CREATE_DATASET(ART_FILES(ARTRES)%FMT, ! FORMAT FICHIER RESULTAT
     *                      ART_FILES(ARTRES)%LU,  ! LU FICHIER RESULTAT
     *                      TITCAS,     ! TITRE DE L'ETUDE
     *                      MAXVAR,     ! MAX VARIABLES SORTIE
     *                      TEXTE,      ! NOMS VARIABLES SORTIE
     *                      SORLEO)     ! SORTIE OU PAS DES VARIABLES
        ! ECRITURE DU MAILLAGE DANS LE FICHIER SORTIE :
        ! SI ON EST ON PARALLEL, FAUT L'INDIQUER VIA NCSIZE ET NPTIR.
        ! LES AUTRES INFORMATIONS SONT DANS MESH.
        ! EN PLUS : DATE/TEMPS DE DEPART ET LES COORDONNEES DE
        ! L'ORIGINE.
        CALL WRITE_MESH(ART_FILES(ARTRES)%FMT, ! FORMAT FICHIER RESULTAT     
     *                  ART_FILES(ARTRES)%LU,  ! LU FICHIER RESULTAT
     *                  MESH,          ! DESCRIPTEUR MAILLAGE
     *                  1,             ! NOMBRE DE PLAN /NA/
     *                  MARDAT,        ! DATE DEBUT
     *                  MARTIM,        ! HEURE DEBUT
     *                  I_ORIG,J_ORIG) ! COORDONNEES DE L'ORIGINE.
  
      
C-----------------------------------------------------------------------
C
C     INITIALISATION DU BLOC PRIVE
C
      IF(NPRIV.GT.0) CALL OS('X=C     ',PRIVE,PRIVE,PRIVE,0.D0)
C
C=======================================================================
C
      IF(NCSIZE.GT.1) THEN
         NFRLIQ=0
         DO I=1,NPTFR
            NFRLIQ=MAX(NFRLIQ,NUMLIQ%I(I))
         ENDDO
         NFRLIQ=P_IMAX(NFRLIQ)
         WRITE(LU,*) ' '
         IF(LNG.EQ.1) WRITE(LU,*) 'NOMBRE DE FRONTIERES LIQUIDES :',
     *        NFRLIQ
         IF(LNG.EQ.2) WRITE(LU,*) 'NUMBER OF LIQUID BOUNDARIES:',NFRLIQ
      ELSE
         CALL FRONT2(NFRLIQ,NFRSOL,DEBLIQ,FINLIQ,DEBSOL,FINSOL,
     *        LIHBOR%I,LIUBOR%I,
     *        MESH%X%R,MESH%Y%R,MESH%NBOR%I,MESH%KP1BOR%I,
     *        IT1%I,NPOIN,NPTFR,KLOG,LISTIN,NUMLIQ%I,MAXFRO)

      ENDIF
C REPERAGE DES FRONTIERES
C
      
C=======================================================================
C
C CORRECTION EVENTUELLE DES VALEURS DU FOND
C
C EN STANDARD, CORFON NE FAIT RIEN
C
      CALL CORFON 
    

C
C=======================================================================
C
C INITIALISATION DE LA HAUTEUR DE HOULE ALEATOIRE A 0.
C
      IF (ALEMON .OR. ALEMUL) THEN
       CALL OS('X=C     ', HALE , SBID , SBID , 0.D0 )
      ENDIF


C
C CALCULS DES DIFFERENTES PERIODES POUR UN CALCUL EN HOULE ALEATOIRE
C
      IF (ALEMON.OR.ALEMUL) THEN
         CALL PERALE(PALE%R,GAMMA,PERPIC,NPALE,T1%R,NPOIN,PRIVE,
     *               NPRIV,PMIN,PMAX)
         PER = PALE%R(1)
      ENDIF
  
      
      

C CALCULS DES DIFFERENTES DIRECTIONS POUR UN CALCUL EN HOULE ALEATOIRE
C MULTIDIRECTIONNELLE
C
      IF (ALEMUL) THEN
         CALL DIRALE(DALE%R,EXPOS,TETAH,TETMIN,TETMAX,NDALE,
     *               T1%R,NPOIN,PRIVE,NPRIV)
      ENDIF
      
     
      

C
C=======================================================================
C
C DEBUT D'UN CALCUL D'AGITATION
C
C LT INDIQUE LE NUMERO DU CALCUL COURANT (ON COMMENCE A 0 AFIN QUE LE
C PREMIER CALCUL SOIT TOUJOURS ENREGISTRE)
C
      LT = 0
C
C INITIALISATION DE QB, T01, T02 ET TM A 0 EN TOUT DEBUT DE CALCUL
C
      CALL OS('X=C     ', QB , SBID , SBID , 0.D0 )
      CALL OS('X=C     ', T01 , SBID , SBID , 0.D0 )
      CALL OS('X=C     ', T02 , SBID , SBID , 0.D0 )
      CALL OS('X=C     ', TM , SBID , SBID , 0.D0 )
C
C
C INITIALISATION DES CONTRAINTES DE RADIATION ET DES
C FORCAGES
C
      CALL OS('X=C     ', FX , SBID , SBID , 0.D0 )
      CALL OS('X=C     ', FY , SBID , SBID , 0.D0 )
      CALL OS('X=C     ', SXX , SBID , SBID , 0.D0 )
      CALL OS('X=C     ', SXY , SBID , SBID , 0.D0 )
      CALL OS('X=C     ', SYY , SBID , SBID , 0.D0 )
      CALL OS('X=C     ', MCOS , SBID , SBID , 0.D0 )
      CALL OS('X=C     ', MSIN , SBID , SBID , 0.D0 )
C
C POUR UN CALCUL EN HOULE ALEATOIRE, LPER ET LDIR INDIQUENT LA PERIODE
C ET LA DIRECTION CALCULEES
C
      LPER = 1
      LDIR = 1
C
100   CONTINUE
C
      IF (BALAYE) THEN
         CALL ENTART(1,PER,LT,LPER,NPERBA,ALEMON,ALEMUL,BALAYE)
      ELSE
         CALL ENTART(1,PER,LT,LPER,NPALE,ALEMON,ALEMUL,BALAYE)
      ENDIF
C
C=======================================================================
C
C : 2                  INITIALISATIONS
C
C=======================================================================
C
C INITIALISATION DES GRANDEURS PHYSIQUES
C     
     

      CALL CONDIH 
C
C=======================================================================
C
C : 3                  CONDITIONS AUX LIMITES
C
C=======================================================================
C
C APPEL DU SOUS-PROGRAMME UTILISATEUR
C


!     
      IF (NCSIZE .GT. 1) THEN 
         CALL BUILD_GLOBAL_BOUND(MESH%KNOLG%I,MESH%NPOIN,NPOIN_TOT,
     c        MESH%NPTFR,NPTFR_TOT,
     c        X,Y,K%R,C%R,CG%R,LIHBOR%I,XT,
     c        YT,KT,CTT,CGT,LIHBORT,MESH%NBOR%I,NBOR_TOT)
      ELSE
           
         DO I=1,NPOIN
            XT(I)=X(I)
            YT(I)=Y(I)
         END DO 
         DO I=1,NPTFR
            NBOR_TOT(I)=MESH%NBOR%I(I)
            LIHBORT(I)=LIHBOR%I(I)
         END DO 
      END IF
     
      CALL BORH
C
C CONSTRUCTION DES MASQUAGES POUR LES CONDITIONS AUX LIMITES
C
     
C EN HOULE ALEATOIRE MULTIDIRECTIONNELLE, LES DIRECTIONS DE PROPAGATION
C AUX LIMITES ONT ETE CALCULEES DANS DALE.
C
200   IF (ALEMUL) THEN
         CALL OS('X=C     ', TETAB ,SBID,SBID, DALE%R(LDIR) )
         CALL ENTART(2,DALE%R(LDIR),LT,LDIR,NDALE,ALEMON,ALEMUL,BALAYE)
      ENDIF
C
C ON CALCULE LES CONDITIONS AUX LIMITES SUR LE POTENTIEL A PARTIR DE CE
C QUE L'UTILISATEUR A DONNE.
C     
c      IF (LT .EQ. 0) THEN

      CALL MASQUE_ARTEMIS
    
  
      CALL PHBOR
c      END IF
C
C=======================================================================
C
C : 4                  RESOLUTION DE L'EQUATION DE BERKHOFF
C
C=======================================================================
C
      CALL BERKHO (LT)
C     
C=======================================================================
C
C : 5.1        CALCUL DES VITESSES, COTE DE LA SURFACE LIBRE, 
C              HAUTEUR ET PHASE DE LA HOULE
C
C=======================================================================
C
      CALL CALRES 
C
      IF (ALEMON .OR. ALEMUL) THEN
C
C        CALCUL CUMULATIF DES MOMENTS m1, m2, et mT1 MIS
C        JUSQU'AU DERNIER CALCUL DANS LES VARIABLES
C        T01, T02, ET TM
C 
         CALL CALCMN
C
      ENDIF
C
C
C=======================================================================
C
C : 5.2        CALCUL DES CONTRAINTES DE RADIATION ET
C              DES FORCES MOTRICES EN HOULE REGULIERE.
C
C=======================================================================
C
      IF (.NOT.ALEMON .AND. .NOT.ALEMUL) THEN
C
       IF (LISHOU) THEN
         CALL DISMOY
     *   (NPOIN,NELEM,MESH%X%R,MESH%Y%R,MESH%IKLE%I,K%R,LISHHO)
       ELSE
         LISHHO = 0
       ENDIF
C
         CALL RADIA1 (LISHHO)
C
      ELSE
         LISHHO = 0
      ENDIF
    
C=======================================================================
C
C : 6   APPEL D'UN SOUS-PROGRAMME UTILISATEUR POUR DES IMPRESSIONS,
C       AU CALCUL DE SOLUTIONS ANALYTIQUES...(NE FAIT RIEN EN STANDARD)
C
C=======================================================================
C
      CALL UTIMP
     *(PHIR%R,PHII%R,C%R,CG%R,K%R,MESH%X%R,MESH%Y%R,ZF%R,H%R,
     * HHO%R,U0%R,V0%R,PHAS%R,S%R,T1%R,T2%R,T3%R,T4%R,INCI%R,
     * GRAV,PER,OMEGA,MESH%IKLE%I,MESH%NBOR%I,MESH%KP1BOR%I,
     * NELEM,NELMAX,IELM,IELMB,NPTFR,NPOIN,PRIVE)
C
C=======================================================================
C
C : 7                  IMPRESSION DES RESULTATS
C
C=======================================================================
C
C
C POUR UN CALCUL EN HOULE ALEATOIRE, LES SORTIES NE SONT FAITES QUE POUR
C LA PERIODE DE PIC
C
      IF (.NOT.ALEMON .AND. .NOT.ALEMUL) THEN
C
C=======================================================================
C
C     CONVERSION DE INCI EN DEGRES
C
C=======================================================================
C
         CALL OS('X=CX    ', INCI , SBID , SBID , RADDEG )
C
C FICHIER RUBENS
C
         CALL BIEF_DESIMP(ART_FILES(ARTRES)%FMT,VARSOR,
     *            HIST,0,NPOIN,ART_FILES(ARTRES)%LU,'STD',PER,0,
     *            LISPRD,LEOPRD,
     *            SORLEO,SORIMP,MAXVAR,TEXTE,0,0)
C
C=======================================================================
C
C              COMPARAISON AVEC UN FICHIER DE REFERENCE
C
C=======================================================================
C
C     LE SOUS-PROGRAMME VALIDA DE LA BIBLIOTHEQUE EST STANDARD.
C     IL PEUT ETRE MODIFIE POUR CHAQUE CAS PARTICULIER.
C     L'APPEL DOIT ETRE LAISSE DANS LA BOUCLE EN TEMPS.
C
         IF(VALID) THEN
           CALL BIEF_VALIDA(TB,TEXTE,
     *                      ART_FILES(ARTREF)%LU,ART_FILES(ARTREF)%FMT,
     *                      VARSOR,TEXTE,
     *                      ART_FILES(ARTRES)%LU,ART_FILES(ARTRES)%FMT,
     *                      MAXVAR,NPOIN,LT,LT,ALIRE)
         ENDIF
C
      ENDIF
C
C=======================================================================
C
C : 8                  PASSAGE A LA PERIODE SUIVANTE
C
C=======================================================================
C
C CAS D'UN BALAYAGE D'UN INTERVALLE DE PERIODE
C
      IF (BALAYE) THEN
         LT   = LT  + 1
         LPER = LPER + 1
         PER  = PER + PERPAS
        
         IF (PER.LE.PERFIN) GOTO 100
      ENDIF
C
C
C=======================================================================
C
C CAS D'UN CALCUL EN HOULE ALEATOIRE
C
C=======================================================================
C
      IF (ALEMON .OR. ALEMUL) THEN
C
         LT  = LT  + 1
C
         IF (LT.LT.NPALE*NDALE) THEN
C
C           ON REACTUALISE L'ENERGIE DE LA HOULE ALEATOIRE
            CALL OS('X=X+CYZ ',HALE,HHO,HHO,1.D0/DBLE(NPALE*NDALE))
C
C           CALCUL SUIVANT EN DIRECTION
            LDIR = LDIR + 1
            IF (LDIR.LE.NDALE) GOTO 200
C
C           CALCUL SUIVANT EN PERIODE
            LDIR = 1
            LPER = LPER + 1
            PER = PALE%R(LPER)
            GOTO 100
C
         ELSE
C
C           DERNIER CALCUL : ON DETERMINE LES VALEURS DES
C           PERIODES MOYENNES (T01 ET T02), ET DE LA DIRECTION
C           MOYENNE (INCI)
C
            CALL CALCTM
C
C           DETERMINATION DES QUANTITES MOYENNES K, C ET CG
C
            CALL CALRE2
C
C           PRISE EN COMPTE DE LA DERNIERE
C           HAUTEUR DE HOULE ALEATOIRE
C
            CALL OS('X=X+CYZ ',HALE,HHO,HHO,1.D0/DBLE(NPALE*NDALE))
            CALL OS('X=SQR(Y)', HALE , HALE , SBID , BID )
            CALL OS('X=CX    ',QB,SBID,SBID,1.D0/DBLE(NPALE*NDALE))
C        
C=======================================================================
C
C           CALCUL DES CONTRAINTES DE RADIATION
C           ET DES FORCES MOTRICES EN HOULE ALEATOIRE
C
C=======================================================================
C
            CALL RADIA2 (LISHHO)
C
C=======================================================================
C
C          CONVERSION DE INCI EN DEGRES
C
C=======================================================================
C
            CALL OS('X=CX    ', INCI , SBID , SBID , RADDEG )
C
C=======================================================================
C
C           FICHIER RUBENS 
C
C=======================================================================
C
            CALL BIEF_DESIMP(ART_FILES(ARTRES)%FMT,VARSOR,
     *            HIST,0,NPOIN,ART_FILES(ARTRES)%LU,'STD',PERPIC,0,
     *            LISPRD,LEOPRD,
     *            SORLEO,SORIMP,MAXVAR,TEXTE,0,0)
C
C=======================================================================
C
C              COMPARAISON AVEC UN FICHIER DE REFERENCE
C
C=======================================================================
C
C     LE SOUS-PROGRAMME VALIDA DE LA BIBLIOTHEQUE EST STANDARD.
C     IL PEUT ETRE MODIFIE POUR CHAQUE CAS PARTICULIER.
C     L'APPEL DOIT ETRE LAISSE DANS LA BOUCLE EN TEMPS.
C
           
            IF(VALID) THEN
              CALL BIEF_VALIDA(TB,TEXTE,
     *                       ART_FILES(ARTREF)%LU,ART_FILES(ARTREF)%FMT,
     *                       VARSOR,TEXTE,
     *                       ART_FILES(ARTRES)%LU,ART_FILES(ARTRES)%FMT,
     *                       MAXVAR,NPOIN,LT,LT,ALIRE)
            ENDIF
C
         ENDIF
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
