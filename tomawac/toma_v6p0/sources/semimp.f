C                       *****************
                        SUBROUTINE SEMIMP
C                       *****************
C
     *( F     , XK    , FREQ  , DFREQ , DEPTH , VENTX , VENTY , X     ,
     *  Y     , NVEB  , NVEF  , NBOR  , NPTFR , DDC   , TV1   , TV2   ,
     *  NP    ,
     *  XRELV , YRELV , U1    , V1    , U2    , V2    , TETA  , SINTET,
     *  COSTET, INDIC , TAILF , RAISF , GRAVIT, CFROT1, CMOUT1, CMOUT2,
     *  TPROP , DTSI  , ROAIR , ROEAU , XKAPPA, BETAM , DECAL , CDRAG ,
     *  ALPHA , ZVENT , NF    , NPLAN , NPOIN2, IANGNL, COEFNL, F1    ,
     *  NSITS , SMOUT , SFROT , SVENT , STRIF , VENT  , VENSTA, VX_CTE,
     *  VY_CTE, SBREK , ALFABJ,
     *  GAMBJ1, GAMBJ2, IQBBJ , IHMBJ , IFRBJ , BORETG, GAMATG, IWHTG ,
     *  IFRTG , ALFARO, GAMARO, GAM2RO, IDISRO, IEXPRO, IFRRO , BETAIH,
     *  EM2SIH, IFRIH , COEFHS, XDTBRK, NDTBRK, STRIA , ALFLTA, RFMLTA,
     *  KSPB  , BDISPB, BDSSPB, PROINF, DF_LIM, LIMIT , CIMPLI, 
     *  NOMVEB, NOMVEF, BINVEN, NBD   , QINDI , TAUWAV,
     *  USOLD , TWOLD , Z0OLD , TSTOT , TSDER , TOLD  , TNEW  , VARIAN,
     *  FMOY  , XKMOY , USNEW , Z0NEW , TWNEW , TAUX1 , TAUX2 , TAUX3 ,
     *  TAUX4 , TAUX5 , TAUX6 , TAUX7 , TRA01 , BETA)
C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  EDF/DER - LABORATOIRE NATIONAL D'HYDRAULIQUE  -  CHATOU (FRANCE)   C
C                                                                     C
C  CODE TOMAWAC DE MODELISATION DES ETATS DE MER EN ELEMENTS FINIS    C
C                       ---- VERSION 5.0 ----                         C
C                                                                     C
C  SEMIMP :  RESOLUTION DE L'ETAPE D'INTEGRATION DES TERMES-SOURCES   C
C  ********  A L'AIDE D'UN SCHEMA A DEGRE D'IMPLICITATION VARIABLE.   C
C                                                                     C
C   - CREE POUR VERSION 1.0  LE 26/03/95 PAR M. BENOIT                C
C   - MOD. POUR VERSION 1.2  LE 07/11/96 PAR M. BENOIT                C
C   - MOD. POUR VERSION 5.0  LE 25/08/00                              C
C
C   - MOD. POUR VERSION 5.9 LE 16/12/2008 PAR JMH
C                           BETA AJOUTE EN ARGUMENT ET PASSE
C                           EN DERNIER ARGUMENT DE QBREK1,2,3,4
C                           A LA PLACE DE TAUX1
C
C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  ARGUMENTS D'APPEL :                                                C
C  *******************                                                C
C  +-------------+----+--------------------------------------------+  C
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !  C
C  +-------------+----+--------------------------------------------+  C
C  ! F(-,-,-)    !<-->! SPECTRE DIRECTIONNEL DE VARIANCE           !  C
C  ! XK(-,-)     ! -->! TABLEAU DES NOMBRES D'ONDE                 !  C
C  ! FREQ(-)     ! -->! TABLEAU DES FREQUENCES DE DISCRETISATION   !  C
C  ! DFREQ(-)    ! -->! TABLEAU DES PAS DE FREQUENCE               !  C
C  ! DEPTH(-)    ! -->! TABLEAU DES PROFONDEURS (METRES)           !  C
C  ! VENTX(-)    !<-->! TABLEAU DE VENT (COMP. OUEST-EST)          !  C
C  ! VENTY(-)    !<-->! TABLEAU DE VENT (COMP. SUD-NORD)           !  C
C  ! X(-)        ! -->! TABLEAU DES ABSCISSES DES POINTS MAILLAGE  !  C
C  ! Y(-)        ! -->! TABLEAU DES ABSCISSES DES POINTS MAILLAGE  !  C
C  ! NVEB-F      ! -->! NUMERO DU FICHIER DES VENTS EN LECTURE     !  C
C  ! NBOR(-)     ! -->! NUMEROS GLOBAUX DES POINTS FRONTIERES      !  C
C  ! NPTFR       ! -->! NOMBRE DE POINTS FRONTIERES                !  C
C  ! DDC         ! -->! DATE DE DEBUT DU CALCUL                    !  C
C  ! TV1         !<-->! DATE DU CHAMP DE CHAMP 1                   !  C
C  ! TV2         !<-->! DATE DU CHAMP DE CHAMP 2                   !  C
C  ! NP          ! -->! NOMBRE DE POINTS DU CHAMP DE VENT          !  C
C  ! XRELV(-)    ! -->! TABLEAU DES ABSCISSES DES POINTS RELEVES   !  C
C  ! YRELV(-)    ! -->! TABLEAU DES ORDONNEES DES POINTS RELEVES   !  C
C  ! U1(-)       !<-->! TABLEAU DES COMP. OUEST-EST DE VENT 1      !  C
C  ! V1(-)       !<-->! TABLEAU DES COMP. SUD-NORD  DE VENT 1      !  C
C  ! U2(-)       !<-->! TABLEAU DES COMP. OUEST-EST DE VENT 2      !  C
C  ! V2(-)       !<-->! TABLEAU DES COMP. SUD-NORD  DE VENT 2      !  C
C  ! TETA(-)     ! -->! VECTEUR DES DIRECTIONS DE DISCRETISATION   !  C
C  ! SINTET(-)   ! -->! VECTEUR DES   SINUS DES DIRECTIONS         !  C
C  ! COSTET(-)   ! -->! VECTEUR DES COSINUS DES DIRECTIONS         !  C
C  ! INDIC       ! -->! TYPE DE FORMAT DE LECTURE                  !  C
C  ! TAILF       ! -->! FACTEUR DE QUEUE DU SPECTRE                !  C
C  ! RAISF       ! -->! RAISON FREQUENTIELLE POUR DISCRETISATION   !  C
C  ! GRAVIT      ! -->! ACCELERATION DE LA PESANTEUR               !  C
C  ! CFROT1      ! -->! CONSTANTE POUR LE TERME DE FROTTEMENT      !  C
C  ! CMOUT1      ! -->! CONSTANTE 1 POUR LE TERME DE MOUTONNEMENT  !  C
C  ! CMOUT2      ! -->! CONSTANTE 2 POUR LE TERME DE MOUTONNEMENT  !  C
C  ! TPROP       ! -->! DATE DE FIN DE L'ETAPE D'INTEGRATION       !  C
C  ! DTSI        ! -->! PAS DE TEMPS D'INTEGRATION (SECONDES)      !  C
C  ! ROAIR       ! -->! MASSE VOLUMIQUE DE L AIR                   !  C
C  ! ROEAU       ! -->! MASSE VOLUMIQUE DE L EAU                   !  C
C  ! XKAPPA      ! -->! CONSTANTE DE VON KARMAN                    !  C
C  ! BETAM       ! -->! CONSTANTE BETAMAX DE LA FORMULE DU VENT    !  C
C  ! DECAL       ! -->! CONSTANTE DE DECALAGE DE CROISSANCE VENT   !  C
C  ! CDRAG       ! -->! COEFFICIENT DE TRAINEE                     !  C
C  ! ALPHA       ! -->! CONSTANTE DE LA LOI DE CHARNOCK            !  C
C  ! ZVENT       ! -->! COTE A LAQUELLE EST MESURE LE VENT (M)     !  C
C  ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !  C
C  ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !  C
C  ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL       !  C
C  ! IANGNL(-,-) ! -->! INDICES ANGULAIRES POUR QUADRUPLETS DIA    !  C
C  ! COEFNL(-)   ! -->! COEFFICIENTS POUR QUADRUPLETS DIA          !  C
C  ! F1          ! -->! PREMIERE FREQUENCE DE DISCRETISATION       !  C
C  ! NSITS       ! -->! NOMBRE DE PAS DE TEMPS D'INTEGRATION       !  C
C  ! SMOUT       ! -->! INDICATEUR DE TYPE DE TERME MOUTONNEMENT   !  C
C  ! SFROT       ! -->! INDICATEUR DE TYPE DE TERME FROTTEMENT     !  C
C  ! SVENT       ! -->! INDICATEUR DE TYPE DE TERME INPUT PAR VENT !  C
C  ! STRIF       ! -->! INDICATEUR DE TYPE DE TERME INTERACTIONS   !  C
C  ! VENT        ! -->! INDICATEUR DE PRISE EN COMPTE DE VENT      !  C
C  ! SBREK       ! -->! INDICATEUR DE TYPE DE TERME DEFERLEMENT    !  C
C  ! ALFABJ      ! -->! MODELE DEFERLEMENT BJ : CONSTANTE ALPHA    !  C
C  ! GAMBJ1      ! -->! MODELE DEFERLEMENT BJ : CONSTANTE GAMMA1   !  C
C  ! GAMBJ2      ! -->! MODELE DEFERLEMENT BJ : CONSTANTE GAMMA2   !  C
C  ! IQBBJ       ! -->! MODELE DEFERLEMENT BJ : MODE CALCUL DE QB  !  C
C  ! IHMBJ       ! -->! MODELE DEFERLEMENT BJ : MODE CALCUL DE HM  !  C
C  ! IFRBJ       ! -->! MODELE DEFERLEMENT BJ : MODE CALCUL DE FREQ!  C
C  ! BORETG      ! -->! MODELE DEFERLEMENT TG : CONSTANTE B        !  C
C  ! GAMATG      ! -->! MODELE DEFERLEMENT TG : CONSTANTE GAMMA    !  C
C  ! IWHTG       ! -->! MODELE DEFERLEMENT TG : MODE CALCUL DE W(H)!  C
C  ! IFRTG       ! -->! MODELE DEFERLEMENT TG : MODE CALCUL DE FREQ!  C
C  ! ALFARO      ! -->! MODELE DEFERLEMENT RO : CONSTANTE ALPHA    !  C
C  ! GAMARO      ! -->! MODELE DEFERLEMENT RO : CONSTANTE GAMMA    !  C
C  ! GAM2RO      ! -->! MODELE DEFERLEMENT RO : CONSTANTE GAMMA2   !  C
C  ! IDISRO      ! -->! MODELE DEFERLEMENT RO : DISTRIBUTION HOULE !  C
C  ! IEXPRO      ! -->! MODELE DEFERLEMENT RO : EXPOSANT N         !  C
C  ! IFRRO       ! -->! MODELE DEFERLEMENT RO : MODE CALCUL DE FREQ!  C
C  ! BETAIH      ! -->! MODELE DEFERLEMENT IH : CONSTANTE BETA     !  C
C  ! EM2SIH      ! -->! MODELE DEFERLEMENT IH : CONSTANTE M2*      !  C
C  ! IFRIH       ! -->! MODELE DEFERLEMENT IH : MODE CALCUL DE FREQ!  C
C  ! COEFHS      ! -->! COEFFICIENT LIMITATEUR DE LA HAUTEUR HS    !  C
C  ! XDTBRK      ! -->! PAS DE TEMPS POUR LE DEFERLEMENT           !  C
C  ! NDTBRK      ! -->! NOMBRE DE SOUS-PAS DE TEMPS DE DEFERLEMENT !  C
C  ! PROINF      ! -->! INDICATEUR DE PROFONDEUR INFINIE           !  C
C  ! NOMVEB-F    ! -->! NOM DU FICHIER DE VENT EN ENTREE           !  C
C  ! BINVEN      ! -->! BINAIRE DU FICHIER DE VENT EN ENTREE       !  C
C  ! TAUWAV(-)   !<-->! TABLEAU DES CONTRAINTES DUES A LA HOULE    !  C
C  ! USOLD(-)    !<-->! VITESSE DE FROTTEMENT A L'INSTANT N        !  C
C  ! TWOLD(-)    !<-->! DIRECTION DU VENT A L'INSTANT N            !  C
C  ! Z0OLD(-)    !<-->! LONGUEUR DE RUGOSITE L'INSTANT N           !  C
C  ! TSTOT(-,-,-)!<-- ! CONTRIBUTION TERME SOURCE - PARTIE TOTALE  !  C
C  ! TSDER(-,-,-)!<-- ! CONTRIBUTION TERME SOURCE - PARTIE DERIVEE !  C
C  ! TOLD(-,-)   !<-->! TABLEAU DE TRAVAIL 1 (DIMENSION NPOIN3)    !  C
C  ! TNEW(-,-)   !<-->! TABLEAU DE TRAVAIL 2 (DIMENSION NPOIN3)    !  C
C  ! VARIAN(-)   !<-- ! TABLEAU DES VARIANCES TOTALES              !  C
C  ! FMOY(-)     !<-- ! TABLEAU DES FREQUENCES MOYENNES            !  C
C  ! XKMOY(-)    !<-- ! TABLEAU DES NOMBRES D'ONDE MOYENS          !  C
C  ! USNEW(-)    !<-->! VITESSE DE FROTTEMENT A L'INSTANT N+1      !  C
C  ! Z0NEW(-)    !<-->! LONGUEUR DE RUGOSITE L'INSTANT N+1         !  C
C  ! TWNEW(-)    !<-->! DIRECTION DU VENT A L'INSTANT N+1          !  C
C  ! TAUX1(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  ! TAUX2(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  ! TAUX3(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  ! TAUX4(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  ! TAUX5(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  ! TAUX6(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  ! TAUX7(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  +-------------+----+--------------------------------------------+  C
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !  C
C  +---------------------------------------------------------------+  C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  APPELS :    - PROGRAMME(S) APPELANT  :  WAC                        C
C  ********    - PROGRAMME(S) APPELE(S) :  USTAR1, TOTNRJ, ANAVEN,    C
C                                          FREMOY, KMOYEN, OV    ,    C
C                                          QWIND1, STRESS, QNLIN1,    C
C                                          QMOUT1, QFROT1, NOUDON,    C
C                                          QWIND2, USTAR2, QBREK1,    C
C                                          QBREK2, QBREK3, QBREK4,    C
C                                          FPREAD, FREM01, FREM02,    C
C                                          FPEPIC                     C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  REMARQUES :                                                        C
C  ***********                                                        C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/ LNG,LU
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER          NPOIN2, NPLAN , NF    , NSITS , NPTFR , NVEB  ,
     *                 NVEF  , LIMIT ,
     *                 SMOUT , SFROT , SVENT , STRIF , SBREK , INDIC ,
     *                 IQBBJ , IHMBJ , IFRBJ , IWHTG , IFRTG , IFRRO ,
     *                 IEXPRO, IFRIH , NDTBRK, NP    , IDISRO, STRIA ,
     *                 NBOR(NPTFR)   , IANGNL(NPLAN,8)
      INTEGER          NBD   , QINDI(NBD)
      DOUBLE PRECISION TAILF , CFROT1, GRAVIT, RAISF , DTSI  , TPROP ,
     *                 CMOUT1, CMOUT2, DDC   , TV1   , TV2   , ZVENT ,
     *                 ROAIR , ROEAU , XKAPPA, BETAM , DECAL , CDRAG ,
     *                 ALPHA , GAMBJ1, GAMBJ2, ALFABJ, BORETG, GAMATG,
     *                 COEFHS, VX_CTE, VY_CTE, CIMPLI,
     *                 GAMARO, ALFARO, GAM2RO, EM2SIH, BETAIH, XDTBRK,
     *                 ALFLTA, RFMLTA, KSPB  , BDISPB, BDSSPB, F1
      DOUBLE PRECISION  DEPTH(NPOIN2), USNEW(NPOIN2) , USOLD(NPOIN2) ,
     *                 VARIAN(NPOIN2),  FMOY(NPOIN2) , XKMOY(NPOIN2) ,
     *                  TWOLD(NPOIN2), TWNEW(NPOIN2) , Z0OLD(NPOIN2) ,
     *                  Z0NEW(NPOIN2), VENTX(NPOIN2) , VENTY(NPOIN2) ,
     *                     U1(NPOIN2),    U2(NPOIN2) ,    V1(NPOIN2) ,
     *                     V2(NPOIN2),     X(NPOIN2) ,     Y(NPOIN2) ,
     *                 TAUWAV(NPOIN2), TAUX1(NPOIN2) , TAUX2(NPOIN2) ,
     *                  TAUX3(NPOIN2), TAUX4(NPOIN2) , TAUX5(NPOIN2) ,
     *                  TAUX6(NPOIN2), TAUX7(NPOIN2) , COEFNL(16)    ,
     *                    TETA(NPLAN), SINTET(NPLAN) , COSTET(NPLAN) ,
     *                    F(NPOIN2,NPLAN,NF),   XK(NPOIN2,NF)        ,
     *                    DF_LIM(NPOIN2,NF) ,
     *                 TSDER(NPOIN2,NPLAN,NF),TSTOT(NPOIN2,NPLAN,NF) ,
     *                 FREQ(NF), DFREQ(NF), XRELV(NP), YRELV(NP)     ,
     *                 TOLD(NPOIN2,NPLAN), TNEW(NPOIN2,NPLAN)
      DOUBLE PRECISION BETA(NPOIN2)
      DOUBLE PRECISION TRA01(NPOIN2,NPLAN)
      CHARACTER*144 NOMVEB, NOMVEF
      CHARACTER*3 BINVEN
      LOGICAL  PROINF, VENT , VENSTA
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER          ISITS , NPOIN3, NPOIN4, IFF   , IP    , JP    ,
     *                 IFCAR , MF1   , MF2   , MFMAX , IDT
      DOUBLE PRECISION AUX1  , AUX2  , AUX3  , AUX4  , COEF  , DFMAX ,
     *                 DEUPI , FM1   , FM2   , TDEB  , TFIN  , VITVEN,
     *                 VITMIN, HM0   , HM0MAX, DTN   , SUM   , AUXI  ,
     *                 USMIN
      CHARACTER*7      CHDON
C
C
      NPOIN3=NPOIN2*NPLAN
      NPOIN4=NPOIN3*NF
      DEUPI=2.D0*3.141592654D0
      VITMIN=1.D-3
C
C
C     -----------------------------------------------------------------
C     ECRETAGE DU SPECTRE EN FONCTION DE LA BATHYMETRIE.
C     -----------------------------------------------------------------
      IF (.NOT.PROINF) THEN
C
C       0.1 CALCUL DE LA VARIANCE TOTALE DU SPECTRE.
C       --------------------------------------------
        CALL TOTNRJ
     *( VARIAN, F     , FREQ  , DFREQ , TAILF , NF    , NPLAN , NPOIN2)
C
C       0.2 CALCUL DU COEFFICIENT DE CORRECTION SUR LE SPECTRE.
C       -------------------------------------------------------
C
        DO IP=1,NPOIN2
          HM0MAX=COEFHS*DEPTH(IP)
          HM0 =MAX(4.D0*DSQRT(VARIAN(IP)),1.D-20)
          TAUX1(IP)=MIN((HM0MAX/HM0)**2,1.D0)
        ENDDO
C
C       0.3 CORRECTION DU SPECTRE.
C       --------------------------
C
        DO IFF=1,NF
          DO JP=1,NPLAN
            DO IP=1,NPOIN2
              F(IP,JP,IFF)=F(IP,JP,IFF)*TAUX1(IP)
            ENDDO
          ENDDO
        ENDDO
C
      ENDIF
C
C     ----------------------------------------------------------------
C     S'IL Y A DU VENT ET QU'IL EST STATIONNAIRE, ON DUPLIQUE LES
C     CONDITIONS DE DEBUT DE PAS DE TEMPS POUR LA FIN DU PAS DE TEMPS.
C     (CAR LES TABLEAUX TWNEW, USNEW ET Z0NEW SONT DES TABLEAUX
C     DE TRAVAIL UTILISES PAR DUMP2D ENTRE 2 APPELS A SEMIMP).
C     ----------------------------------------------------------------
      IF (VENT.AND.VENSTA) THEN
        DO IP=1,NPOIN2
          TWNEW(IP)=TWOLD(IP)
        ENDDO
        IF (SVENT.EQ.2) THEN
          DO IP=1,NPOIN2
            USNEW(IP)=USOLD(IP)
            Z0NEW(IP)=Z0OLD(IP)
          ENDDO
        ENDIF
      ENDIF
C
C
C     -----------------------------------------------------------------
C     DEBUT DE BOUCLE PRINCIPALE SUR LE NOMBRE (NSITS) DE PAS DE TEMPS
C     D INTEGRATION DES TERMES-SOURCES PAR PAS DE TEMPS DE PROPAGATION.
C     -----------------------------------------------------------------
      DO 100 ISITS=1,NSITS
C
C
C       1. AFFECTATION DES DATES DE DEBUT ET FIN DE PAS DE TEMPS.
C       =========================================================
        TDEB=TPROP-DBLE(NSITS-ISITS+1)*DTSI
        TFIN=TDEB+DTSI
C
C
C       2. MISE A JOUR (EVENTUELLE) DES TABLEAUX DE VENT.
C       =================================================
        IF (VENT.AND..NOT.VENSTA) THEN
C
C         2.1 MISE A JOUR DU CHAMP DE VENT POUR LA DATE TFIN.
C         ---------------------------------------------------
          CHDON='VENT   '
          IF (NOMVEB(1:1).NE.' ') THEN
            CALL NOUDON
     *( VENTX , VENTY , X     , Y     , NPOIN2, NVEB  , BINVEN, NBOR  ,
     *  NPTFR , TFIN  , DDC   , TV1   , TV2   , NP    , XRELV , YRELV ,
     *  TOLD  , TNEW  , TRA01 , U1    , V1    , U2    , V2    ,
     *  INDIC , CHDON , 2 ) 
          ELSEIF (NOMVEF(1:1).NE.' ') THEN
            CALL NOUDON
     *( VENTX , VENTY , X     , Y     , NPOIN2, NVEF  , BINVEN, NBOR  ,
     *  NPTFR , TFIN  , DDC   , TV1   , TV2   , NP    , XRELV , YRELV ,
     *  TOLD  , TNEW  , TRA01 , U1    , V1    , U2    , V2    ,
     *  INDIC , CHDON , 2 )
          ELSE
            CALL ANAVEN
     *( VENTX , VENTY , X     , Y     , NPOIN2, TFIN  , DDC   , VX_CTE,
     *  VY_CTE)
          ENDIF
C
C         2.2 CALCUL DE LA DIRECTION DU VENT.
C         -----------------------------------
C
          DO IP=1,NPOIN2
            VITVEN=SQRT(VENTX(IP)**2+VENTY(IP)**2)
            IF (VITVEN.GT.VITMIN) THEN
              TWNEW(IP)=ATAN2(VENTX(IP),VENTY(IP))
            ELSE
              TWNEW(IP)=0.D0
            ENDIF
          ENDDO
C
C         2.3 CALCUL DES VITESSES DE FROTTEMENT ET LONGUEURS RUGOSITE.
C         ------------------------------------------------------------
          IF (SVENT.EQ.2) CALL USTAR2
     *( USNEW , VENTX , VENTY , NPOIN2)
C
        ENDIF
C
        IF (VENT) THEN
          IF (SVENT.EQ.1) CALL USTAR1
     *( USNEW , Z0NEW , TAUWAV, VENTX , VENTY , CDRAG , ALPHA , XKAPPA,
     *  ZVENT , GRAVIT, NPOIN2)
        ENDIF
C
C
C       3. CALCUL DES GRANDEURS MOYENNES DU SPECTRE DIRECTIONNEL.
C       =========================================================
C
C       3.1 CALCUL DE LA VARIANCE TOTALE DU SPECTRE.
C       --------------------------------------------
        CALL TOTNRJ
     *( VARIAN, F     , FREQ  , DFREQ , TAILF , NF    , NPLAN , NPOIN2)
C
C       3.2 CALCUL DE LA FREQUENCE MOYENNE DU SPECTRE.
C       ----------------------------------------------
        CALL FREMOY
     *( FMOY  , F     , FREQ  , DFREQ , TAILF , NF    , NPLAN , NPOIN2,
     *  TAUX1 , TAUX2 )
C
C       3.3 CALCUL DU NOMBRE D'ONDE MOYEN DU SPECTRE.
C       ---------------------------------------------
        CALL KMOYEN
     *( XKMOY , XK    , F     , FREQ  , DFREQ , TAILF , NF    , NPLAN ,
     *  NPOIN2, TAUX1 , TAUX2 , TAUX3 )
C
C
C       4. CALCUL DES CONTRIBUTIONS DES TERMES-SOURCES DE GENERATION,
C          DE MOUTONNEMENT ET D'INTERACTIONS ENTRE QUADRUPLETS.
C       =============================================================
C
C       4.1 INITIALISATIONS DES TABLEAUX DES TERMES-SOURCES.
C       ----------------------------------------------------
        DO IFF=1,NF
          DO JP=1,NPLAN
            DO IP=1,NPOIN2
              TSTOT(IP,JP,IFF)=0.0D0
              TSDER(IP,JP,IFF)=0.0D0
            ENDDO
          ENDDO
        ENDDO
C
C       4.2 GENERATION PAR LE VENT.
C       ---------------------------
        IF (VENT) THEN
          IF (SVENT.EQ.1) THEN
            CALL QWIND1
     *( TSTOT , TSDER , F     , XK    , FREQ  , USOLD , USNEW , TWOLD ,
     *  TWNEW , Z0OLD , Z0NEW , TETA  , ROAIR , ROEAU , BETAM , XKAPPA,
     *  DECAL , GRAVIT, NF    , NPLAN , NPOIN2, CIMPLI, TOLD  , TNEW  ,
     *  TAUX1 , TAUX2 , TAUX3 , TAUX4 , TAUX5 , TAUX6 , TAUX7 )
            CALL STRESS
     *( TAUWAV, TSTOT , F     , USNEW , TWNEW , Z0NEW , FREQ  , DFREQ ,
     *  TETA  , SINTET, COSTET, ROAIR , ROEAU , XKAPPA, BETAM , DECAL ,
     *  GRAVIT, NPOIN2, NPLAN , NF    , TAUX1 , TAUX2 , TAUX3 )
          ELSE
            CALL QWIND2
     *( TSTOT , TSDER , F     , XK    , FREQ  , USOLD , USNEW , TWOLD ,
     *  TWNEW , TETA  , ROAIR , ROEAU , GRAVIT, NF    , NPLAN , NPOIN2,
     *  CIMPLI, TAUX1 , TAUX2 , TAUX3 , TAUX4 , TAUX5 )
          ENDIF
        ELSE
          DO IP=1,NPOIN2
            USNEW(IP)=0.D0
          ENDDO
        ENDIF
C
C       4.3 INTERACTIONS NON-LINEAIRES ENTRE QUADRUPLETS FREQUENTIELS.
C       --------------------------------------------------------------
        IF (STRIF.EQ.1) CALL QNLIN1
     *( TSTOT , TSDER , IANGNL, COEFNL, NF    , NPLAN , F1    , RAISF ,
     *  TAILF , PROINF, NPOIN2, F     , DEPTH , XKMOY , TAUX1 , TAUX2 ,
     *  TAUX3 , TAUX4 , TAUX5 , TAUX6 )
C
C
C       4.4 DISSIPATION PAR MOUTONNEMENT (WHITE-CAPPING).
C       -------------------------------------------------
        IF (SMOUT.EQ.1) CALL QMOUT1
     *( TSTOT , TSDER , F     , XK    , VARIAN, FREQ  , FMOY  , XKMOY ,
     *  PROINF, CMOUT1, CMOUT2, GRAVIT, NF    , NPLAN , NPOIN2, TAUX1 ,
     *  TAUX2 )
C
C       4.5 DISSIPATION PAR FROTTEMENT SUR LE FOND.
C       -------------------------------------------
        IF ((SFROT.EQ.1).AND.(.NOT.PROINF)) CALL QFROT1
     *( TSTOT , TSDER , F     , XK    , DEPTH , CFROT1, GRAVIT, NF    ,
     *  NPLAN , NPOIN2, TAUX1 )

C.......4.6 CALCUL DU LIMITEUR DE CROISSANCE
C       ------------------------------------
C.......Pas de limiteur (en fait on met une valeur tres elevee).
        IF (LIMIT.EQ.0) THEN
          AUXI=1.D99
          DO IFF=1,NF
            DO IP=1,NPOIN2
              DF_LIM(IP,IFF)=AUXI
            ENDDO
          ENDDO
C
C.......Limiteur de la version cycle 4 de WAM.
        ELSEIF (LIMIT.EQ.1) THEN
C          COEF=6.4D-7*GRAVIT**2*DTSI/1200.D0	
          COEF=0.62D-4*DTSI/1200.D0
          DO IFF=1,NF
            AUXI=COEF/FREQ(IFF)**5
            DO IP=1,NPOIN2
              DF_LIM(IP,IFF)=AUXI
            ENDDO
          ENDDO
C
C.......Limiteur de Hersbach et Janssen (1999) sans la partie Uetoile.
        ELSEIF (LIMIT.EQ.2) THEN
          COEF=3.0D-7*GRAVIT*FREQ(NF)*DTSI
          DO IFF=1,NF
            AUXI=COEF/FREQ(IFF)**4
            USMIN=GRAVIT*5.6D-3/FREQ(IFF)
            DO IP=1,NPOIN2
              DF_LIM(IP,IFF)=AUXI*MAX(USNEW(IP),USMIN)
            ENDDO
          ENDDO    
        ENDIF  
C
C
C       5. MISE A JOUR SPECTRE - PRISE EN COMPTE DES TERMES-SOURCES DE
C          GENERATION, DE MOUTONNEMENT ET D'INTERACTIONS QUADRUPLETS.
C       ==============================================================
C
        DO IFF=1,NF
          DO JP=1,NPLAN
            DO IP=1,NPOIN2
              DFMAX=DF_LIM(IP,IFF)
              AUX1=MAX( 1.D0-DTSI*TSDER(IP,JP,IFF)*CIMPLI , 1.D0 )
              AUX2=DTSI*TSTOT(IP,JP,IFF)/AUX1
              AUX3=MIN( ABS(AUX2) , DFMAX )
              AUX4=SIGN(AUX3,AUX2)
              F(IP,JP,IFF)=MAX( F(IP,JP,IFF)+AUX4 , 0.D0 )
            ENDDO
          ENDDO
        ENDDO
C
C
C       6. TRAITEMENT SPECIAL POUR LA PARTIE HAUTES-FREQUENCES.
C       =======================================================
C
C       6.1 CALCUL DE LA FREQUENCE MOYENNE DU SPECTRE.
C       ----------------------------------------------
        CALL FREMOY
     *( FMOY  , F     , FREQ  , DFREQ , TAILF , NF    , NPLAN , NPOIN2,
     *  TAUX1 , TAUX2 )
C
        AUX1=GRAVIT/(7.D0*DEUPI*FREQ(1))
        AUX2=2.5D0/FREQ(1)
        AUX3=1.D0/LOG10(RAISF)
C
        DO IP=1,NPOIN2
C
C       6.2 CALCUL DE LA DERNIERE FREQUENCE DE LA PARTIE DISCRETISEE.
C           CETTE FREQUENCE EST DEFINIE COMME LE MAXIMUM DE
C           FM1=4.*FPM   ET  FM2=2.5*FMOY. SON INDICE EST NOTE MFMAX.
C       -------------------------------------------------------------
          FM1 =AUX1/(USNEW(IP)+1.D-90)
          FM2 =AUX2*FMOY(IP)
          MF1=INT(AUX3*LOG10(FM1)+1.D0)
          MF2=INT(AUX3*LOG10(FM2)+1.D0)
          MFMAX=MIN( MAX(MF1,MF2) , NF )
C
C       6.3 MODIFICATION DE LA PARTIE HAUTES FREQUENCES.
C           UNE DECROISSANCE EN F**(-TAILF) EST IMPOSEE
C           AU-DELA DE FREQ(MFMAX).  (TAILF=5. DANS WAM CYCLE 4)
C       --------------------------------------------------------
          DO IFF=MFMAX+1,NF
            AUX4=(FREQ(MFMAX)/FREQ(IFF))**TAILF
            DO JP=1,NPLAN
              F(IP,JP,IFF)=AUX4*F(IP,JP,MFMAX)
            ENDDO
          ENDDO
        ENDDO
C
C       7. PRISE EN COMPTE DU TERME-SOURCE DE DEFERLEMENT
C       =================================================
C
        IF (((SBREK.GT.0).OR.(STRIA.GT.0)).AND.(.NOT.PROINF)) THEN 
C
C         7.1 CALCUL D'UNE FREQUENCE REPRESENTATIVE.
C         ------------------------------------------
          IF ((SBREK.GT.0).AND.(SBREK.LT.5)) THEN
            IF (SBREK.EQ.1) IFCAR = IFRBJ
            IF (SBREK.EQ.2) IFCAR = IFRTG
            IF (SBREK.GE.3) IFCAR = IFRRO
            IF (SBREK.GE.4) IFCAR = IFRIH
C
            GOTO (751,752,753,754,755,756), IFCAR
            IF(LNG.EQ.1) THEN
              WRITE(LU,*) 'FREQUENCE DE HOULE NON PREVUE......IFCAR=',
     *                     IFCAR
            ELSE
              WRITE(LU,*) 'WAVE FREQUENCY NOT EXPECTED......IFCAR=',
     *                     IFCAR
            ENDIF
            GOTO 759
C
C           FREQUENCE MOYENNE FMOY.
C           - - - - - - - - - - - -
  751       CONTINUE
            DO IP=1,NPOIN2
              TAUX3(IP)=FMOY(IP)
            ENDDO
            GOTO 759
C
C           FREQUENCE MOYENNE F01.
C           - - - - - - - - - - -
  752       CONTINUE
            CALL FREM01
     *( TAUX3 , F     , FREQ  , DFREQ , TAILF , NF    , NPLAN , NPOIN2,
     *  TAUX1 , TAUX2 )
            GOTO 759
C
C           FREQUENCE MOYENNE F02.
C           - - - - - - - - - - -
  753       CONTINUE
            CALL FREM02
     *( TAUX3 , F     , FREQ  , DFREQ , TAILF , NF    , NPLAN , NPOIN2,
     *  TAUX1 , TAUX2 )
            GOTO 759
C
C           FREQUENCE DE PIC (FREQUENCE DISCRETE PORTANT MAX VARIANCE).
C           - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  754       CONTINUE
            CALL FREPIC
     *( TAUX3 , F     , FREQ  , NF    , NPLAN , NPOIN2, TAUX1 , TAUX2 )
            GOTO 759
C
C           FREQUENCE DE PIC (METHODE DE READ AVEC EXPOSANT 5).
C           - - - - - - - - - - - - - - - - - - - - - - - - - - 
  755       CONTINUE
            CALL FPREAD
     *( TAUX3 , F     , FREQ  , DFREQ , NF    , NPLAN , NPOIN2, 5.D0  ,
     *  TAILF , TAUX1 , TAUX2 )
            GOTO 759
C
C           FREQUENCE DE PIC (METHODE DE READ AVEC EXPOSANT 8).
C           - - - - - - - - - - - - - - - - - - - - - - - - - - 
  756       CONTINUE
            CALL FPREAD
     *( TAUX3 , F     , FREQ  , DFREQ , NF    , NPLAN , NPOIN2, 8.D0  ,
     *  TAILF , TAUX1 , TAUX2 )
C
  759       CONTINUE
C
        ENDIF
C
C.........BOUCLE SUR LES SOUS-PAS DE TEMPS POUR LE DEFERLEMENT.
C         = = = = = = = = = = = = = = = = = = = = = = = = = = =
          SUM=(XDTBRK**NDTBRK-1.D0)/(XDTBRK-1.D0)
          DTN=DTSI/SUM
C
          DO 782 IDT=1,NDTBRK

C         7.2 INITIALISATIONS DES TABLEAUX DES TERMES-SOURCES.
C         ----------------------------------------------------
          DO IFF=1,NF
            DO JP=1,NPLAN
              DO IP=1,NPOIN2
                TSTOT(IP,JP,IFF)=0.0D0
              ENDDO
            ENDDO
          ENDDO
C
C         7.3 CALCUL DE LA VARIANCE TOTALE DU SPECTRE.
C         --------------------------------------------
          CALL TOTNRJ
     *( VARIAN, F     , FREQ  , DFREQ , TAILF , NF    , NPLAN , NPOIN2)
C
C         7.4 CALCUL CONTRIBUTION DU DEFERLEMENT
C         --------------------------------------
          GOTO (761,762,763,764), SBREK
          IF(SBREK.NE.0) THEN
           IF(LNG.EQ.1) THEN
           WRITE(LU,*) 'TYPE DE DEFERLEMENT NON IMPLANTE...SBREK=',SBREK
           WRITE(LU,*) 'PAS DE PRISE EN COMPTE DU DEFERLEMENT'
           ELSE
           WRITE(LU,*) 'BREAKING FORMULATION NOT PROGRAMMED...SBREK=',
     *                  SBREK
           WRITE(LU,*) 'NO CONSIDERATION OF THE DEPTH-INDUCED BREAKING'
           ENDIF
          ENDIF
          GOTO 769
C
C         7.4.1 MODELE DE DEFERLEMENT DE BATTJES ET JANSSEN (1978).
C         - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  761     CONTINUE
          CALL QBREK1
     *( TSTOT , TSDER , F     , TAUX3 , VARIAN, DEPTH , ALFABJ, GAMBJ1,
     *  GAMBJ2, IQBBJ , IHMBJ , NF    , NPLAN , NPOIN2, BETA )
          GOTO 769
C
C         7.4.2 MODELE DE DEFERLEMENT DE THORNTON ET GUZA (1983).
C         - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  762     CONTINUE
          CALL QBREK2
     *( TSTOT , TSDER , F     , TAUX3 , VARIAN, DEPTH , BORETG, GAMATG,
     *  IWHTG , NF    , NPLAN , NPOIN2, BETA )
          GOTO 769
C
C         7.4.3 MODELE DE DEFERLEMENT ROELVINK (1993).
C         - - - - - - - - - - - - - - - - - - - - - -
  763     CONTINUE
          CALL QBREK3
     *( TSTOT , TSDER , F     , TAUX3 , VARIAN, DEPTH , ALFARO, GAMARO,
     *  GAM2RO, IEXPRO, IDISRO, NF    , NPLAN , NPOIN2, BETA )
          GOTO 769
C
C         7.4.4 MODELE DE DEFERLEMENT IZUMIYA ET HORIKAWA (1984).
C         - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  764     CONTINUE
          CALL QBREK4
     *( TSTOT , TSDER , F     , TAUX3 , VARIAN, DEPTH , BETAIH, EM2SIH,
     *  GRAVIT, NF    , NPLAN , NPOIN2, BETA )
          GOTO 769
C
  769     CONTINUE
C
C       7.5 INTERACTIONS NON-LINEAIRES ENTRE TRIPLETS FREQUENTIELS.
C       -----------------------------------------------------------
          IF (STRIA.EQ.1) THEN
            CALL FREMOY
     *( FMOY  , F     , FREQ  , DFREQ , TAILF , NF    , NPLAN , NPOIN2,
     *  TAUX1 , TAUX2 )
            CALL QTRIA1
     *( F     , XK    , FREQ  , DEPTH , RAISF , GRAVIT, ALFLTA, RFMLTA,
     *  NF    , NPLAN , NPOIN2, TSTOT , TSDER , VARIAN, FMOY  )
C
        ELSEIF (STRIA.EQ.2) THEN
            CALL QTRIA2
     *( F     , XK    , FREQ  , DFREQ , DEPTH , TETA  , SINTET, COSTET ,
     *  KSPB  , BDISPB, BDSSPB, RAISF , GRAVIT, NF    , NPLAN , NPOIN2 ,
     *  NBD   , QINDI , TSTOT , TSDER )
        ENDIF
C
C         7.5 MISE A JOUR SPECTRE - PRISE EN COMPTE DU TERME-SOURCE
C             DE DEFERLEMENT (SCHEMA EXPLICITE D'EULER)
C         ---------------------------------------------------------
C
        DO IFF=1,NF
          DO JP=1,NPLAN
            DO IP=1,NPOIN2
              F(IP,JP,IFF)=MAX(F(IP,JP,IFF)+DTN*TSTOT(IP,JP,IFF),0.D0)
            ENDDO
          ENDDO
        ENDDO
C
        DTN=DTN*XDTBRK
C
  782   CONTINUE
C
        ENDIF
C
C
C       8. PASSAGE DES TABLEAUX NEW EN OLD POUR PROCHAIN PAS DE TEMPS.
C       ==============================================================
        IF (VENT) THEN
          DO IP=1,NPOIN2
            USOLD(IP)=USNEW(IP)
            Z0OLD(IP)=Z0NEW(IP)
            TWOLD(IP)=TWNEW(IP)
          ENDDO
        ENDIF
C
C
  100 CONTINUE
C     -----------------------------------------------------------------
C     FIN DE LA BOUCLE PRINCIPALE SUR LE NOMBRE (NSITS) DE PAS DE TEMPS
C     D INTEGRATION DES TERMES-SOURCES PAR PAS DE TEMPS DE PROPAGATION.
C     -----------------------------------------------------------------
C
      RETURN
      END
