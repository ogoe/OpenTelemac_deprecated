C                       *****************
                        SUBROUTINE DUMP2D
C                       *****************
C
     *( LT , DEUPI , XF1 , NP1 )
C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  EDF/DER - LABORATOIRE NATIONAL D'HYDRAULIQUE  -  CHATOU (FRANCE)   C
C                                                                     C
C  CODE TOMAWAC DE MODELISATION DES ETATS DE MER EN ELEMENTS FINIS    C
C                       ---- VERSION 1.2 ----                         C
C                                                                     C
C  DUMP2D :  IMPRESSION DES VARIABLES SYNTHETIQUES DE HOULE, DE VENT, C
C  ********  DE COURANT, DE BATHYMETRIE, ... AUX NOEUDS DU MAILLAGE   C
C            SPATIAL 2D (FORMAT SELAFIN BINAIRE)                      C
C                                                                     C
C   - CREE POUR VERSION 1.0  LE 01/02/95 PAR F. MARCOS                C
C   - MOD. POUR VERSION 1.2  LE 04/07/96 PAR M. BENOIT                C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  ARGUMENTS D'APPEL :                                                C
C  *******************                                                C
C  +-------------+----+--------------------------------------------+  C
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !  C
C  +-------------+----+--------------------------------------------+  C
C  ! F(-,-,-)    ! -->! SPECTRE DIRECTIONNEL DE VARIANCE           !  C
C  ! AT          ! -->! DATE COURANTE DU CALCUL                    !  C
C  ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL 2D    !  C
C  ! NELEM2      ! -->! NOMBRE D'ELEMENTS DU MAILLAGE SPATIAL 2D   !  C
C  ! IKLE2(-,-)  ! -->! CORRESPONDANCE NUMEROTATION LOCALE-GLOBALE !  C
C  ! NBOR(-)     ! -->! NUMEROTATION DES POINTS FRONTIERE          !  C
C  ! NPTFR       ! -->! NOMBRE DE  POINTS FRONTIERE                !  C
C  ! ITR01       !<-->! TABLEAU DE TRAVAIL DE DIMENSION 3*NELEM2   !  C
C  ! SORG2D      ! -->! INDICATEUR DE SORTIE DES VARIABLES 2D      !  C
C  ! X(-)        ! -->! TABLEAU DES ABSCISSES DES POINTS MAILLAGE  !  C
C  ! Y(-)        ! -->! TABLEAU DES ABSCISSES DES POINTS MAILLAGE  !  C
C  ! DEPTH(-)    ! -->! TABLEAU DES PROFONDEURS (METRES)           !  C
C  ! UX(-)       ! -->! TABLEAU DES COMP. OUEST-EST DU COURANT     !  C
C  ! UY(-)       ! -->! TABLEAU DES COMP. SUD-NORD  DU COURANT     !  C
C  ! VENTX(-)    !<-->! TABLEAU DE VENT (COMP. OUEST-EST)          !  C
C  ! VENTY(-)    !<-->! TABLEAU DE VENT (COMP. SUD-NORD)           !  C
C  ! COURAN      ! -->! INDICATEUR DE PRESENCE D'UN COURAN         !  C
C  ! PROINF      ! -->! INDICATEUR DE PROFONDEUR INFINIE           !  C
C  ! VENT        ! -->! INDICATEUR DE PRESENCE D'UN VENT           !  C
C  ! NR2D        ! -->! NUM. DU FICHIER DE SORTIE DES RESULTATS 2D !  C
C  ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !  C
C  ! TAILF       ! -->! FACTEUR DE QUEUE DU SPECTRE                !  C
C  ! BINR2D      ! -->! TYPE DE BINAIRE DU FICHIER DES RESULTATS 2D!  C
C  ! TITCAS      ! -->! TITRE DU CAS DE CALCUL                     !  C
C  ! PRIVE(-,-,-)! -->! TABLEAU UTILISATEUR                        !  C
C  ! NPRIV       ! -->! DIMENSION DU TABLEAU UTILISATEUR           !  C
C  ! TRA02(-,-)  !<-->! TABLEAU DE TRAVAIL                         !  C
C  ! TRA03(-,-)  !<-->! TABLEAU DE TRAVAIL                         !  C
C  ! TRA04(-,-)  !<-->! TABLEAU DE TRAVAIL CONTENANT U*,Z0,...     !  C
C  ! DEBRES      ! -->! INDICATEUR DE PREMIERE DATE A SAUVER       !  C
C  ! XK(-,-)     ! -->! TABLEAU DES NOMBRES D'ONDE                 !  C
C  ! CG(-,-)     ! -->! TABLEAU DES VITESSES DE GROUPE             !  C
C  ! W1(-,-)     !<-->! TABLEAU DE TRAVAIL                         !  C
C  ! MESH(-)     ! -->!                                            !  C
C  ! XMESH(-)    ! -->! STRUCTURE MAILLAGE                         !  C
C  ! COSTET(-)   ! -->! VECTEUR DES COSINUS DES DIRECTIONS         !  C
C  ! SINTET(-)   ! -->! VECTEUR DES SINUS   DES DIRECTIONS         !  C
C  ! ROEAU       ! -->! MASSE VOLUMIQUE DE L'EAU                   !  C
C  ! GRAVIT      ! -->! ACCELERATION DE LA PESANTEUR               !  C
C  ! DEUPI       ! -->! 2.PI                                       !  C
C  ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !  C
C  ! FREQ(-)     ! -->! TABLEAU DES FREQUENCES DE DISCRETISATION   !  C
C  ! DFREQ(-)    ! -->! TABLEAU DES PAS DE FREQUENCES              !  C
C  ! IELM2       ! -->!                                            !  C
C  ! T1(-)       !<-->! TABLEAU DE TRAVAIL STRUCTURE               !  C
C  ! T2(-)       !<-->! TABLEAU DE TRAVAIL STRUCTURE               !  C
C  ! T3(-)       !<-->! TABLEAU DE TRAVAIL STRUCTURE               !  C
C  ! T4(-)       !<-->! TABLEAU DE TRAVAIL STRUCTURE               !  C
C  +-------------+----+--------------------------------------------+  C
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !  C
C  +---------------------------------------------------------------+  C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  APPELS :    - PROGRAMME(S) APPELANT  :  TOMWAC                     C
C  ********    - PROGRAMME(S) APPELE(S) :  ECRIT                      C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C
      USE DECLARATIONS_TOMAWAC
C
      IMPLICIT NONE
C
      INTEGER          LT    , I1    , I2    , NP1   , IP
      DOUBLE PRECISION DEUPI , RADDEG, PIS2  , U10   , FMIN  , FMAX
      DOUBLE PRECISION XF1(NP1)
C
C
      RADDEG=57.29577951D0
      PIS2=1.570796327D0
      FMIN=FREQ(1)
      FMAX=FREQ(NF)
C
C
C=====C
C     C CALCUL DES VARIABLES SELECTIONNEES.
C=====C====================================
C L'ordre de calcul des variables ne correspond pas à celui des sorties
C graphiques afin de limiter le nombre de tableaux de travail.
C
C     -------------------------------CONTRAINTES DE RADIATION
      I1=NPOIN3*NF+1
      I2=I1+NPOIN2*NF
      IF (.NOT.PROINF) THEN
        IF ( SORLEO(11).OR.SORLEO(12).OR.SORLEO(13).OR.
     *       SORLEO(14).OR.SORLEO(15) ) CALL RADIAT
     *( STRA51%R, STRA52%R, STRA53%R, STRA54%R, STRA55%R,
     *  SXK%R   , XF1        , SCG%R   , SDEPTH%R,
     *  TRA02(I1:I2),STRA36%R, STRA37%R, STRA38%R, STRA39%R,
     *  NPOIN2)
      ENDIF
C     -------------------------------ETALEMENT DIRECTIONNEL
      IF (SORLEO(4)) THEN
        CALL SPREAD
     *( STRA31%R, XF1        , SCOSTE%R, SSINTE%R, NPLAN ,
     *  SFR%R   , SDFR%R  , NF         , NPOIN2     , TAILF      ,
     *  STRA34%R, STRA35%R, STRA36%R, STRA37%R, STRA38%R,
     *  STRA39%R)
      ENDIF
C     -------------------------------DIRECTION MOYENNE
      IF (SORLEO(3)) THEN
       CALL TETMOY
     *( STRA32%R, XF1   , SCOSTE%R, SSINTE%R, NPLAN , FREQ  ,
     *  SDFR%R  , NF    , NPOIN2     , TAILF      , STRA36%R   ,
     *  STRA37%R, STRA38%R, STRA39%R )
       IF (TRIGO) THEN
         DO IP=1,NPOIN2
           TRA32(IP)=(PIS2-TRA32(IP))*RADDEG
         ENDDO
       ELSE
         DO IP=1,NPOIN2
           TRA32(IP)=TRA32(IP)*RADDEG
         ENDDO
       ENDIF
      ENDIF
C     -------------------------------FREQUENCE MOYENNE FMOY
      IF (SORLEO(18).OR.SORLEO(28)) THEN
        CALL FREMOY
     *( STRA33%R, XF1   , SFR%R   , SDFR%R  , TAILF , NF  ,
     *  NPLAN      , NPOIN2, STRA38%R, STRA39%R)
        IF (SORLEO(28)) THEN
          DO IP=1,NPOIN2
            TRA61(IP)=1.D0/MIN(MAX(TRA33(IP),FMIN),FMAX)
          ENDDO
        ENDIF
      ENDIF
C     -------------------------------FREQUENCE MOYENNE FM01
      IF (SORLEO(19).OR.SORLEO(29)) THEN
        CALL FREM01
     *( STRA34%R, XF1   , SFR%R   , SDFR%R  , TAILF , NF  ,
     *  NPLAN      , NPOIN2, STRA38%R, STRA39%R)
        IF (SORLEO(29)) THEN
          DO IP=1,NPOIN2
            TRA62(IP)=1.D0/MIN(MAX(TRA34(IP),FMIN),FMAX)
          ENDDO
        ENDIF
      ENDIF
C     -------------------------------FREQUENCE MOYENNE FM02
      IF (SORLEO(20).OR.SORLEO(30)) THEN
        CALL FREM02
     *( STRA35%R, XF1   , SFR%R   , SDFR%R  , TAILF , NF  ,
     *  NPLAN      , NPOIN2, STRA38%R, STRA39%R)
        IF (SORLEO(30)) THEN
          DO IP=1,NPOIN2
            TRA63(IP)=1.D0/MIN(MAX(TRA35(IP),FMIN),FMAX)
          ENDDO
        ENDIF
      ENDIF
C     -------------------------------FREQUENCE PIC DISCRETE
      IF (SORLEO(21).OR.SORLEO(31)) THEN
        CALL FREPIC
     *( STRA36%R, XF1   , SFR%R , NF   , NPLAN , NPOIN2,
     *  STRA38%R, STRA39%R      )
        IF (SORLEO(31)) THEN
          DO IP=1,NPOIN2
            TRA64(IP)=1.D0/MIN(MAX(TRA36(IP),FMIN),FMAX)
          ENDDO
        ENDIF
      ENDIF
C     -------------------------------FREQUENCE PIC READ ORDRE 5
      IF (SORLEO(22).OR.SORLEO(32)) THEN
        CALL FPREAD
     *( STRA56%R, XF1   , SFR%R, SDFR%R  , NF   , NPLAN ,
     *  NPOIN2     , 5.D0  , TAILF   , STRA38%R, STRA39%R  )
        IF (SORLEO(32)) THEN
          DO IP=1,NPOIN2
            TRA65(IP)=1.D0/MIN(MAX(TRA56(IP),FMIN),FMAX)
          ENDDO
        ENDIF
      ENDIF
C     -------------------------------FREQUENCE PIC READ ORDRE 8
      IF (SORLEO(23).OR.SORLEO(33)) THEN
        CALL FPREAD
     *( STRA57%R, XF1   , SFR%R , SDFR%R  , NF   , NPLAN ,
     *  NPOIN2     , 8.D0  , TAILF    , STRA38%R, STRA39%R  )
        IF (SORLEO(33)) THEN
          DO IP=1,NPOIN2
            TRA66(IP)=1.D0/MIN(MAX(TRA57(IP),FMIN),FMAX)
          ENDDO
        ENDIF
      ENDIF
C
      IF (VENT) THEN
C       -------------------------------COEFFICIENT DE TRAINEE
        IF (SORLEO(25)) THEN
          DO IP=1,NPOIN2
            U10=UV(IP)**2+VV(IP)**2
            IF (U10.GT.1.D-6) THEN
              TRA58(IP)=TRA42(IP)**2/U10
            ELSE
              TRA58(IP)=0.D0
            ENDIF
          ENDDO
        ENDIF
      ENDIF
C       -------------------------------VITESSE AU FOND
      IF (.NOT.PROINF) THEN
        IF (SORLEO(16)) THEN
          CALL VITFON
     *( STRA59%R, XF1   , SXK%R , SDEPTH%R, SDFR%R , NF   ,
     *  NPOIN2     , NPLAN , GRAVIT   , STRA39%R)
        ENDIF
      ENDIF
C     -------------------------------VARIANCE
      IF (SORLEO(1).OR.SORLEO(2)) THEN
        CALL TOTNRJ
     *( STRA37%R , XF1   , SFR%R  , SDFR%R , TAILF ,
     *  NF  , NPLAN , NPOIN2)
C     -------------------------------HAUTEUR SIGNIFICATIVE
        IF (SORLEO(2)) THEN
          DO IP=1,NPOIN2
            TRA38(IP)=4.D0*SQRT(TRA37(IP))
          ENDDO
        ENDIF
      ENDIF
C     -------------------------------PUISSANCE LINEIQUE
      IF (SORLEO(34)) THEN
        CALL WPOWER
     *( STRA60%R, XF1   , SFR%R  , SDFR%R , SCG%R  , TAILF , NF   ,
     *  NPLAN , NPOIN2, ROEAU , GRAVIT)
      ENDIF
C
      RETURN
      END
